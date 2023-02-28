//! dna sketching
//! 
//! 
//! 

use std::time::{SystemTime};
use cpu_time::{ProcessTime,ThreadTime};


use std::path::{Path, PathBuf};

// for multithreading
use crossbeam_channel::*;
use serde::{de::DeserializeOwned, Serialize};

use std::fmt::{Debug};

// our crate
use hnsw_rs::prelude::*;


use kmerutils::base::{kmergenerator::*, Kmer32bit, Kmer64bit, CompressedKmerT};
use kmerutils::sketching::seqsketchjaccard::{SeqSketcherT, ProbHash3aSketch, SuperHashSketch};

use crate::utils::{idsketch::*, reloadhnsw};
use crate::utils::files::{process_dir,process_dir_parallel, process_files_group, ProcessingState, DataType};

use crate::dna::dnafiles::{process_file_in_one_block, process_buffer_in_one_block, process_file_concat_split};

use crate::utils::parameters::*;


// Sig is basic item of a signature , VecSig is a vector of such items
type VecSig<Sketcher, Kmer>  = Vec< <Sketcher as SeqSketcherT<Kmer>>::Sig>;

fn sketchandstore_dir_compressedkmer<Kmer:CompressedKmerT+KmerBuilder<Kmer>, Sketcher : SeqSketcherT<Kmer> + Send + Sync>(dirpath : &Path, sketcher : Sketcher ,
                    filter_params: &FilterParams, 
                    processing_params : &ProcessingParams, 
                    other_params : &ComputingParams) 

        where  <Sketcher as SeqSketcherT<Kmer>>::Sig : 'static + Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>, 
                DistHamming : Distance< <Sketcher as  SeqSketcherT<Kmer>>::Sig>,
                {
    //
    //
    log::trace!("sketchandstore_dir processing dir {}", dirpath.to_str().unwrap());
    log::info!("sketchandstore_dir {}", dirpath.to_str().unwrap());
    let start_t = SystemTime::now();
    let cpu_start = ProcessTime::now();
    //
    let block_processing = processing_params.get_block_flag();
    // a queue of signature waiting to be inserted , size must be sufficient to benefit from threaded probminhash and insert
    // and not too large to spare memory
    let insertion_block_size = 5000;
    // set parallel_files to false for backaward compatibility, nb_files_by_group possibly needs to be adjusted 
    let nb_files_by_group = 500;
    let mut insertion_queue : Vec<IdSeq>= Vec::with_capacity(insertion_block_size);
    //
    let mut hnsw : Hnsw::< <Sketcher as SeqSketcherT<Kmer>>::Sig, DistHamming>;
    let mut state :ProcessingState;
    //
    if other_params.get_adding_mode() {
        // in this case we must reload
        let dirpath = std::env::current_dir();
        if dirpath.is_err() {
            log::error!("dnasketch::sketchandstore_dir_compressedkmer cannot get current directory");
            std::panic!("dnasketch::sketchandstore_dir_compressedkmer cannot get current directory");
        }
        let dirpath = dirpath.unwrap();
        log::info!("dnasketch::sketchandstore_dir_compressedkmer will reload hnsw data from director {:?}", dirpath);
        let hnsw_opt = reloadhnsw::reload_hnsw(&dirpath, &AnnParameters::default());
        if hnsw_opt.is_none() {
            log::error!("cannot reload hnsw from directory : {:?}", &dirpath);
            std::process::exit(1);
        }
        else { 
            hnsw = hnsw_opt.unwrap();
        }
        let reload_res = ProcessingState::reload_json(&dirpath);
        if reload_res.is_ok() {
            state = reload_res.unwrap();
        }
        else {
            log::error!("dnasketch::cannot reload processing state (file 'processing_state.json' from directory : {:?}", &dirpath);
            std::process::exit(1);           
        }
    } 
    else {
        // creation mode
        let hnsw_params = processing_params.get_hnsw_params();
        hnsw = Hnsw::< <Sketcher as SeqSketcherT<Kmer>>::Sig, DistHamming>::new(hnsw_params.get_max_nb_connection() as usize , hnsw_params.capacity , 16, hnsw_params.get_ef(), DistHamming{});
        state = ProcessingState::new();
    }
    hnsw.set_extend_candidates(true);
    hnsw.set_keeping_pruned(false);
    //
    // Sketcher allocation, we do not need reverse complement hashing as we sketch assembled genomes. (Jianshu Zhao)
    //
    let kmer_hash_fn = | kmer : &Kmer | -> Kmer::Val {
        let mask : Kmer::Val = num::NumCast::from::<u64>((0b1 << 2*kmer.get_nb_base()) - 1).unwrap();
        let hashval = kmer.get_compressed_value() & mask;
        hashval
    };
//    let sketcher_params = processing_params.get_sketching_params();
//    let sketcher = seqsketchjaccard::SeqSketcher::new(sketcher_params.get_kmer_size(), sketcher_params.get_sketch_size());
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(insertion_block_size);
    // launch process_dir in a thread or async
    crossbeam_utils::thread::scope(|scope| {
        // sequence sending, productor thread
        let sender_handle = scope.spawn(move |_|   {
        let nb_sent : usize;
        let start_t_prod = SystemTime::now();
            let res_nb_sent;
            if block_processing {
                log::info!("dnasketch::sketchandstore_dir_compressedkmer : block processing");
                if other_params.get_parallel_io() {
                    let mut nb_sent_parallel;
                    log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir_parallel");
                    let mut path_block = Vec::<PathBuf>::with_capacity(nb_files_by_group);
                    let res_nb_sent_parallel = process_dir_parallel(&mut state, &DataType::DNA,  dirpath, filter_params, 
                                    nb_files_by_group,  &mut path_block, &process_buffer_in_one_block, &send);
                    // we must treat residue in path_block if any
                    match res_nb_sent_parallel {
                        Ok(nb) => { nb_sent_parallel = nb;},
                        _             => {  log::error!("\n some error occurred in process_dir_parallel");
                                            std::panic!("\n some error occurred in process_dir_parallel");
                                        },
                    };
                    if path_block.len() > 0 {
                        log::info!("dnasketch::process_dir_parallel sending residue, size : {}", path_block.len());
                        let seqs = process_files_group(&DataType::DNA, filter_params, &path_block, &process_buffer_in_one_block);
                        for mut seqfile in seqs {
                            for i in 0..seqfile.len() {
                                seqfile[i].rank = state.nb_seq;
                                state.nb_seq += 1;
                            }                
                            state.nb_file += 1;
                            nb_sent_parallel += 1;
                            send.send(seqfile).unwrap();
                        }
                        path_block.clear();
                        if log::log_enabled!(log::Level::Info) && state.nb_file % 1000 == 0 {
                            log::info!("nb file processed : {}, nb sequences processed : {}", state.nb_file, state.nb_seq);
                        }
                        if state.nb_file % 1000 == 0 {
                            println!("nb file processed : {}, nb sequences processed : {}", state.nb_file, state.nb_seq);
                        }
                    }
                    res_nb_sent = Ok(nb_sent_parallel);
                } // end case parallel
                else {
                    log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir serial");
                    res_nb_sent = process_dir(&mut state, &DataType::DNA,  dirpath, filter_params, 
                                    &process_file_in_one_block, &send);
                }
            }
            else {
                log::info!("processing by concat and split");
                res_nb_sent = process_dir(&mut state, &DataType::DNA, dirpath, filter_params, &process_file_concat_split, &send);
            }
            match res_nb_sent {
                Ok(nb_really_sent) => {
                    nb_sent = nb_really_sent;
                    println!("process_dir processed nb sequences : {}", nb_sent);
                }
                Err(_) => {
                    nb_sent = 0;
                    println!("some error occured in process_dir");
                    log::error!("some error occured in process_dir");
                }
            };
            drop(send);
            state.elapsed_t =  start_t_prod.elapsed().unwrap().as_secs() as f32;
            log::info!("sender processed in  system time(s) : {}", state.elapsed_t);
            // dump processing state in the current directory
            let _ = state.dump_json(&Path::new("./"));
            Box::new(nb_sent)
        });
        //
        // sequence reception, consumer thread
        //
        let receptor_handle = scope.spawn(move |_| {
            let sender_cpu = ThreadTime::try_now();
            let mut seqdict : SeqDict;
            if other_params.get_adding_mode() {
                // must reload seqdict
                let mut filepath = PathBuf::new();
                filepath.push("seqdict.json");
                let res_reload = SeqDict::reload_json(&filepath);
                if res_reload.is_err() {
                    let cwd = std::env::current_dir();
                    if cwd.is_ok() {
                        log::info!("current directory : {:?}", cwd.unwrap());
                    }
                    log::error!("cannot reload SeqDict (file 'seq.json' from current directory");
                    std::process::exit(1);   
                }
                else {
                    seqdict = res_reload.unwrap();
                }
            }
            else {
                seqdict =  SeqDict::new(100000);
            }
            // we must read messages, sketch and insert into hnsw
            let mut read_more = true;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = receive.recv();
                match res_receive {
                    Err(RecvError) => { read_more = false;
                        // sketch the content of  insertion_queue if not empty
                        if insertion_queue.len() > 0 {
                            log::debug!("end of reception");
                            let sequencegroup_ref : Vec<&Sequence> = insertion_queue.iter().map(|s| s.get_sequence_dna().unwrap()).collect();
                            log::debug!("end of reception received nb seq : {}", sequencegroup_ref.len());
                            // collect rank
                            let seq_rank :  Vec<usize> = insertion_queue.iter().map(|s| s.get_rank()).collect();
                            // collect Id
                            let mut seq_id :  Vec<ItemDict> = insertion_queue.iter().map(|s| ItemDict::new(Id::new(s.get_path(), s.get_fasta_id()), s.get_seq_len())).collect();
                            seqdict.0.append(&mut seq_id);
                            let signatures = sketcher.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);                            
                            // we have Vec<u64> signatures we must go back to a vector of IdSketch for hnsw insertion
                            let mut data_for_hnsw = Vec::<(&VecSig<Sketcher,Kmer>, usize)>::with_capacity(signatures.len());
                            for i in 0..signatures.len() {
                                data_for_hnsw.push((&signatures[i], seq_rank[i]));
                            }
                            // parallel insertion
                            hnsw.parallel_insert(&data_for_hnsw);
                        }
                    }
                    Ok(mut idsequences) => {
                        // concat the new idsketch in insertion queue.
                        insertion_queue.append(&mut idsequences);
                        // if insertion_queue is beyond threshold size we can go to threaded sketching and threading insertion
                        if insertion_queue.len() > insertion_block_size {
                            let sequencegroup_ref : Vec<&Sequence> = insertion_queue.iter().map(|s| s.get_sequence_dna().unwrap()).collect();
                            log::debug!("received nb seq : {}", sequencegroup_ref.len());
                            let seq_rank :  Vec<usize> = insertion_queue.iter().map(|s| s.get_rank()).collect();
                            // collect Id
                            let mut seq_id :  Vec<ItemDict> = insertion_queue.iter().map(|s| ItemDict::new(Id::new(s.get_path(), s.get_fasta_id()), s.get_seq_len())).collect();
                            seqdict.0.append(&mut seq_id);
                            // computes hash signature
                            log::debug!("calling sketch_probminhash3a_compressedkmer");
                            let signatures = sketcher.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch, inserting unique id, for hnsw insertion
                            let mut data_for_hnsw = Vec::<(&VecSig<Sketcher,Kmer>, usize)>::with_capacity(signatures.len());
                            for i in 0..signatures.len() {
                                data_for_hnsw.push((&signatures[i], seq_rank[i]));
                            }
                            // parallel insertion
                            hnsw.parallel_insert(&data_for_hnsw);
                            // we reset insertion_queue
                            insertion_queue.clear();
                        }
                    }
                }
            }
            //
            // We must dump hnsw to save "database" if not empty
            //
            if  hnsw.get_nb_point() > 0 {
                let hnswdumpname = String::from("hnswdump");
                log::info!("going to dump hnsw");
                let resdump = hnsw.file_dump(&hnswdumpname);
                match resdump {
                    Err(msg) => {
                        println!("dump failed error msg : {}", msg);
                    },
                    _ =>  { println!("dump of hnsw ended");}
                };
                // dump some info on layer structure
                hnsw.dump_layer_info();
                // dumping dictionary
                let resdump = seqdict.dump(String::from("seqdict.json"));
                match resdump {
                    Err(msg) => {
                        println!("seqdict dump failed error msg : {}", msg);
                    },
                    _ =>  { println!("dump of seqdict ended OK");}
                };                
            }
            else {
                log::info!("no dumping hnsw, no data points");
            }
            // and finally dump processing parameters in file name "parameters.json"
            let _ = processing_params.dump_json(&Path::new("./"));
            // get time for io and fasta parsing
            if sender_cpu.is_ok() {
                let cpu_time = sender_cpu.unwrap().try_elapsed();
                if cpu_time.is_ok() {
                    log::info!("sender needed : {}", cpu_time.unwrap().as_secs());
                }
            }
            //
            Box::new(seqdict.0.len())
        }); // end of receptor thread
        // now we must join handles
        let nb_sent = sender_handle.join().unwrap();
        let nb_received = receptor_handle.join().unwrap();
        log::info!("sketchandstore, nb_sent = {}, nb_received = {}", nb_sent, nb_received);
        if nb_sent != nb_received {
            log::warn!("an error occurred  nb msg sent : {}, nb msg received : {}", nb_sent, nb_received);
        }
    }).unwrap();  // end of scope
    // get total time 
    let cpu_time = cpu_start.elapsed().as_secs();
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;

    if log::log_enabled!(log::Level::Info) {
        log::info!("processing of directory  : total (io+hashing+hnsw) cpu time(s) {}", cpu_time);
        log::info!("processing of directory : total (io+hashing+hnsw) elapsed time(s) {}", elapsed_t);
    }
    else {
        println!("process_dir : cpu time(s) {}", cpu_time);
        println!("process_dir : elapsed time(s) {}", elapsed_t);
    }
} // end of sketchandstore_dir_compressedkmer 


// hashing function  used in sketching, we do not need reverse complement hashing as we sketch assembled genomes. (Jianshu Zhao)
#[allow(unused)]
pub (crate) fn kmer_hash<Kmer>( kmer : &Kmer) -> Kmer::Val 
    where Kmer : CompressedKmerT + KmerBuilder<Kmer> {
    //
    let mask : Kmer::Val = num::NumCast::from::<u64>((0b1 << 2*kmer.get_nb_base()) - 1).unwrap();
    let hashval = kmer.get_compressed_value() & mask;
    hashval
}


/// This drives sequence sketching and sending to hnsw
pub fn dna_process_tohnsw(dirpath : &Path, filter_params : &FilterParams, processing_parameters : &ProcessingParams, others_params : &ComputingParams) {
    // dispatch according to kmer_size
    let kmer_size = processing_parameters.get_sketching_params().get_kmer_size();
    //
    let sketchparams = processing_parameters.get_sketching_params();
    match sketchparams.get_algo() {
        SketchAlgo::PROB3A => {
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = ProbHash3aSketch::<Kmer32bit>::new(sketchparams);
                //
                sketchandstore_dir_compressedkmer::<Kmer32bit, ProbHash3aSketch::<Kmer32bit> >(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = ProbHash3aSketch::<Kmer16b32bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, ProbHash3aSketch::<Kmer16b32bit>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = ProbHash3aSketch::<Kmer64bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer64bit, ProbHash3aSketch::<Kmer64bit>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        } // end PPROB
        SketchAlgo::SUPER => {
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = SuperHashSketch::<Kmer32bit>::new(sketchparams);
                //
                sketchandstore_dir_compressedkmer::<Kmer32bit, SuperHashSketch::<Kmer32bit> >(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = SuperHashSketch::<Kmer16b32bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, SuperHashSketch::<Kmer16b32bit>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = SuperHashSketch::<Kmer64bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer64bit, SuperHashSketch::<Kmer64bit>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        }
    }
} // end of dna_process