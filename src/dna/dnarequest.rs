//! gathers dna request stuf
//! 
//! 


use std::path::{Path, PathBuf};
use std::fs::{OpenOptions};
use std::io::{BufWriter};

use std::time::{SystemTime};
use cpu_time::ProcessTime;

use serde::{Serialize, de::DeserializeOwned};
use std::fmt::{Debug};

#[allow(unused)]
use fxhash::*;

use hnsw_rs::prelude::*;

use kmerutils::base::{kmergenerator::*, Kmer32bit, Kmer64bit, CompressedKmerT};
use kmerutils::sketching::seqsketchjaccard::*;

use probminhash::{setsketcher::SetSketchParams};

use crate::utils::*;
use crate::dna::dnafiles::{process_file_by_sequence, process_buffer_by_sequence, process_buffer_in_one_block, process_file_in_one_block};
use crate::{matcher::*, answer::ReqAnswer};



// Sig is basic item of a signature , VecSig is a vector of such items
type VecSig<Sketcher, Kmer>  = Vec< <Sketcher as SeqSketcherT<Kmer>>::Sig>;


// a type to describe msessage to collector task

struct RequestMsg<Sig> {
    // signature and sequence rank in requests
    pub sketch_and_rank : Vec<Sig>,
    // id of sequence
    pub item : ItemDict,
}


impl <Sig> RequestMsg<Sig> {

    pub fn new(sketch_and_rank : Vec<Sig>, item : ItemDict) -> Self {
        RequestMsg{sketch_and_rank, item}
    }

} // end of RequestMsg






fn sketch_and_request_dir_compressedkmer<Kmer:CompressedKmerT + KmerBuilder<Kmer>,Sketcher : SeqSketcherT<Kmer> + Send + Sync>(request_dirpath : &Path, 
                sketcher : Sketcher,
                filter_params: &FilterParams, seqdict : &SeqDict, 
                processing_parameters : &ProcessingParams,
                other_params : &ComputingParams, 
                hnsw : &Hnsw< <Sketcher as SeqSketcherT<Kmer>>::Sig ,DistHamming>, 
                knbn : usize, ef_search : usize) -> Matcher

            where Kmer::Val : num::PrimInt + Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
                  KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>, 
                  DistHamming : Distance< <Sketcher as SeqSketcherT<Kmer>>::Sig > {
    //
    let sketcher_params = processing_parameters.get_sketching_params();
    let block_processing = processing_parameters.get_block_flag();
    //
    log::trace!("sketch_and_request_dir processing dir {}", request_dirpath.to_str().unwrap());
    log::info!("Dna mode sketch_and_request_dir {}", request_dirpath.to_str().unwrap());
    log::info!("sketch_and_request kmer size  {}  sketch size {} ", sketcher_params.get_kmer_size(), sketcher_params.get_sketch_size());
    let out_threshold = 0.99;  // TODO threshold needs a test to get initialized!
    // creating an output file in the current directory
    let outname = "gsearch.answers.txt";
    let outpath = PathBuf::from(outname.clone());
    let outfile = OpenOptions::new().write(true).create(true).truncate(true).open(&outpath);
    if outfile.is_err() {
        log::error!("sketch_and_request_dir_compressedkmer could not open file {:?}", outpath.as_os_str());
        println!("SeqDict dump: could not open file {:?}", outpath.as_os_str());
        return Err("SeqDict Deserializer dump failed").unwrap();
    }
    let mut outfile = BufWriter::new(outfile.unwrap());
    log::info!("dumping request answers in : {}, thresold dist : {} ", outname, out_threshold);
    //
    let start_t = SystemTime::now();
    let cpu_start = ProcessTime::now();
    // a queue of signature request , size must be sufficient to benefit from threaded probminhash and search
    let request_block_size = 500;
    let mut state = ProcessingState::new();
    
    //
    // Sketcher allocation, we do need reverse complement
    //
    let kmer_hash_fn = | kmer : &Kmer | -> Kmer::Val {
        let canonical =  kmer.reverse_complement().min(*kmer);
        let mask : Kmer::Val = num::NumCast::from::<u64>((0b1 << 2*kmer.get_nb_base()) - 1).unwrap();
        let hashval = canonical.get_compressed_value() & mask;
        hashval
    };
    // create something for likelyhood computation
    let mut matcher = Matcher::new(processing_parameters.get_kmer_size(), sketcher_params.get_sketch_size(), seqdict);
    let mut nb_sent = 0;
    let nb_received = 0;
    //
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(request_block_size + 10);
    // to send sketch result to a collector task
    let (request_sender , request_receiver) = 
        crossbeam_channel::bounded::<RequestMsg<<Sketcher as SeqSketcherT<Kmer>>::Sig>>(request_block_size+1);
    //
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let pool_nb_thread = pool.current_num_threads();
    log::info!("nb threads in pool : {:?}", pool_nb_thread);
    // launch process_dir in a thread or async
    pool.scope(|scope| {
        // sequence sending, productor thread
        scope.spawn(|_|   {
            let res_nb_sent;
            match block_processing {
                true => {
                    log::info!("dnarequest::sketchandstore_dir_compressedkmer : block processing");
                    if other_params.get_parallel_io() {
                        let nb_files_by_group = other_params.get_nb_files_par();
                        log::info!("dnarequest::sketchandstore_dir_compressedkmer : calling process_dir_parallel, nb_files in parallel : {}", nb_files_by_group);
                        res_nb_sent = process_dir_parallel(&mut state, &DataType::DNA, request_dirpath, filter_params, 
                                        nb_files_by_group, &process_buffer_in_one_block, &send);
                    } // end case parallel io
                    else {
                        res_nb_sent = process_dir(&mut state, &DataType::DNA, request_dirpath, &filter_params, &process_file_in_one_block, &send);
                    }
                },
                false => {
                    log::info!("request::sketchandstore_dir_compressedkmer : seq by seq processing");
                    if other_params.get_parallel_io() {
                        let nb_files_by_group = other_params.get_nb_files_par();
                        log::info!("dnarequest::sketchandstore_dir_compressedkmer : calling process_dir_parallel, nb_files in parallel : {}", nb_files_by_group);
                        res_nb_sent = process_dir_parallel(&mut state, &DataType::DNA,  &request_dirpath, filter_params, 
                                nb_files_by_group, &process_buffer_by_sequence, &send);
                    }
                    else {
                        log::info!("processing by sequence, sketching whole file globally");
                        res_nb_sent = process_dir(&mut state, &DataType::DNA, &request_dirpath, filter_params,
                                     &process_file_by_sequence, &send);
                    }
                },
            };
            //
            match res_nb_sent {
                Ok(nb_really_sent) => {
                    nb_sent = nb_really_sent;
                    println!("process_dir processed nb sequences : {}", nb_sent);
                    log::info!("sender processed nb sequences {}", nb_sent);
                }
                Err(_) => {
                    println!("some error occured in process_dir");
                }
            };
            state.elapsed_t =  start_t.elapsed().unwrap().as_secs() as f32;
            log::info!("sender processed in  system time(s) : {}", state.elapsed_t);
            drop(send);
        }); // end of request sender 

        // sequence reception, sketching/request consumer thread
        scope.spawn( |scope| {
            // we must read messages, sketch and insert into hnsw
            let mut sketcher_queue : Vec<IdSeq>= Vec::with_capacity(request_block_size);
            let mut nb_sketched = 0;
            let mut read_more = true;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = receive.recv();
                match res_receive {
                    Err(_) => { read_more = false;
                        log::debug!("end of request reception");
                    },
                    Ok(mut idsequences) => {
                        // concat the new idsketch in insertion queue.
                        log::debug!("request reception nb request : {}", idsequences.len());
                        sketcher_queue.append(&mut idsequences);
                    }
                };
                // if request_queue is beyond threshold size we can go to threaded sketching and threading insertion
                if sketcher_queue.len() > request_block_size  || !read_more {
                    match block_processing { 
                        true => {
                            let sequencegroup_ref : Vec<&Sequence> = sketcher_queue.iter().map(|s| s.get_sequence_dna().unwrap()).collect();
                            // collect Id
                            let seq_item : Vec<ItemDict> = sketcher_queue.iter().map(|s| ItemDict::new(Id::new(s.get_path(), s.get_fasta_id()), s.get_seq_len())).collect();
                            // computes hash signature
                            let signatures = sketcher.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
                            // now we must send signatures, rak and seq_items to request collector
                            for i in 0..signatures.len() {
                                let request_msg = RequestMsg::new(signatures[i].clone() , seq_item[i].clone());
                                let _ = request_sender.send(request_msg);
                            }
                            // update state
                            nb_sketched += signatures.len();
                            sketcher_queue.clear();                            
                        },
                        //
                        false => { // means we are in seq by seq mode inside a file. We get one signature by msg (or file)
                            let sequencegroup_ref : Vec<&Sequence> = sketcher_queue.iter().map(|s| s.get_sequence_dna().unwrap()).collect();
                            // collect Id
                            let seq_item : Vec<ItemDict> = sketcher_queue.iter().map(|s| ItemDict::new(Id::new(s.get_path(), s.get_fasta_id()), s.get_seq_len())).collect();
                            // computes hash signature
                            let signatures = sketcher.sketch_compressedkmer_seqs(&sequencegroup_ref, kmer_hash_fn);
                            // now we must send signatures, rak and seq_items to request collector
                            for i in 0..signatures.len() {
                                let request_msg = RequestMsg::new(signatures[i].clone() , seq_item[i].clone());
                                let _ = request_sender.send(request_msg);
                            }
                            // update state
                            nb_sketched += signatures.len();
                            sketcher_queue.clear();  
                        },
                    };
                }
            } // end while 
            //
            log::info!("nb request sketched : {}", nb_sketched);
        }); // end of receptor sketcher thread
        //
        // now we spawn a new thread that is devoted to is devoted to hnsw search
        //
        scope.spawn(|_| {
            let mut request_store = Vec::<VecSig<Sketcher,Kmer>>::with_capacity(3 * request_block_size);
            let mut itemv = Vec::<ItemDict>::with_capacity(3 * request_block_size);
            let mut read_more = true;
            let mut nb_request = 0;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = request_receiver.recv();
                match res_receive {
                    Err(_recv_error) =>  { read_more = false;
                                                    log::debug!("end of collector reception");
                    }
                    Ok(to_insert) => {
                        request_store.push(to_insert.sketch_and_rank);
                        itemv.push(to_insert.item);
                        nb_request += 1;
                        log::debug!("request collector received nb_received : {}", nb_received);
                    }
                }
                if read_more == false || request_store.len() > request_block_size {
                    log::debug!("inserting block in hnsw, nb new points : {:?}", request_store.len());
                    let mut data_for_hnsw = Vec::<VecSig<Sketcher,Kmer> >::with_capacity(request_store.len());
                    for i in 0..request_store.len() {
                        data_for_hnsw.push(request_store[i].clone());
                    }
                    let knn_neighbours = hnsw.parallel_search(&data_for_hnsw, knbn, ef_search);  
                    // construct and dump answers
                    for i in 0..knn_neighbours.len() {
                        let answer = ReqAnswer::new(nb_request+i, itemv[i].clone(), &knn_neighbours[i]);
                        if answer.dump(&seqdict, out_threshold, &mut outfile).is_err() {
                            log::info!("could not dump answer for request id {}", answer.get_request_id().get_id().get_fasta_id());
                        }
                        // store in matcher. remind that each i corresponds to a request
                        let candidates = knn_neighbours[i].iter().map(|n| SequenceMatch::new(seqdict.0[n.d_id].clone(), n.get_distance())).collect();
                        matcher.insert_sequence_match(itemv[i].clone(), candidates);
                    }
                    nb_request += request_store.len();
                    request_store.clear();
                    itemv.clear();
                }
            }
            //
            log::debug!("request collector thread dumping hnsw , received nb_received : {}", nb_received);
        }); // end of collector thread
    });  // end of pool.scope
    //
    log::debug!("sketch_and_request_dir, nb_sent = {}, nb_received = {}", nb_sent, &nb_received);
    if nb_sent != nb_received {
        log::error!("an error occurred  nb msg sent : {}, nb msg received : {}", nb_sent, nb_received);
    }
    log::info!("matcher collected {} answers", matcher.get_nb_sequence_match());
    let cpu_time = cpu_start.elapsed().as_secs();
    log::info!("process_dir : cpu time(s) {}", cpu_time);
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
    log::info!("process_dir : elapsed time(s) {}", elapsed_t);
    //
    matcher
} // end of sketch_and_request_dir_compressedkmer 





// This function returns paired sequence by probminhash and hnsw 
pub fn get_sequence_matcher(request_params : &RequestParams, processing_params : &ProcessingParams,
                     filter_params : &FilterParams, computing_params : &ComputingParams, seqdict : &SeqDict, 
                     ef_search : usize) -> Result<Matcher, String> {
    //
    let sketch_params = processing_params.get_sketching_params();
    let request_dirpath = Path::new(request_params.get_req_dir());
    let database_dirpath = Path::new(request_params.get_hnsw_dir());
    let nbng = request_params.get_nb_answers();
    log::info!("sketch params reloaded kmer size : {}, sketch size {}", sketch_params.get_kmer_size(), sketch_params.get_sketch_size());
    //
    let matcher : Matcher;
    // reload hnsw
    log::info!("\n reloading hnsw from {}", database_dirpath.to_str().unwrap());

    // allocate sketcher and get matcher
    match sketch_params.get_algo() {
        SketchAlgo::PROB3A => {
            match sketch_params.get_kmer_size() {
                1..=14 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <Kmer32bit as CompressedKmerT>::Val>(database_dirpath)?;
                    let sketcher = ProbHash3aSketch::<Kmer32bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, ProbHash3aSketch::<Kmer32bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &computing_params, 
                            &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <Kmer64bit as CompressedKmerT>::Val>(database_dirpath)?;
                    let sketcher = ProbHash3aSketch::<Kmer64bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, ProbHash3aSketch::<Kmer64bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &computing_params,
                            &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <Kmer16b32bit as CompressedKmerT>::Val>(database_dirpath)?;
                    let sketcher = ProbHash3aSketch::<Kmer16b32bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, ProbHash3aSketch::<Kmer16b32bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &computing_params,
                            &hnsw, nbng as usize, ef_search);
                }
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));
                }
            }
        } //end match PROBA3
        SketchAlgo::SUPER => {
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHashSketch<Kmer32bit, f32> as SeqSketcherT<Kmer32bit> >::Sig >(database_dirpath)?;
                    let sketcher = SuperHashSketch::<Kmer32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, SuperHashSketch::<Kmer32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHashSketch<Kmer64bit, f32> as SeqSketcherT<Kmer64bit> >::Sig >(database_dirpath)?;
                    let sketcher = SuperHashSketch::<Kmer64bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, SuperHashSketch::<Kmer64bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHashSketch<Kmer16b32bit, f32> as SeqSketcherT<Kmer16b32bit> >::Sig >(database_dirpath)?;
                    let sketcher = SuperHashSketch::<Kmer16b32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, SuperHashSketch::<Kmer16b32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search); 
                } 
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));                   
                }          
            }
        }  // end match SUPER
        SketchAlgo::HLL => {
            // 16 bits is sufficient up to kmer of size 30!
            let mut hll_params = SetSketchParams::default();
            if hll_params.get_m() < sketch_params.get_sketch_size() as u64 {
                log::warn!("!!!!!!!!!!!!!!!!!!!  need to adjust hll parameters!");
            }
            hll_params.set_m(sketch_params.get_sketch_size());
            //
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <HyperLogLogSketch<Kmer32bit, u16> as SeqSketcherT<Kmer32bit> >::Sig >(database_dirpath)?;
                    let sketcher = HyperLogLogSketch::<Kmer32bit, u16>::new(sketch_params, hll_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, HyperLogLogSketch::<Kmer32bit, u16> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <HyperLogLogSketch<Kmer64bit, u16> as SeqSketcherT<Kmer64bit> >::Sig >(database_dirpath)?;
                    let sketcher = HyperLogLogSketch::<Kmer64bit, u16>::new(sketch_params, hll_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, HyperLogLogSketch::<Kmer64bit, u16> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <HyperLogLogSketch<Kmer16b32bit, u16> as SeqSketcherT<Kmer16b32bit> >::Sig >(database_dirpath)?;
                    let sketcher = HyperLogLogSketch::<Kmer16b32bit, u16>::new(sketch_params, hll_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, HyperLogLogSketch::<Kmer16b32bit, u16> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search); 
                } 
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));                   
                }          
            }
        } // end match HLL
        SketchAlgo::SUPER2 => {
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHash2Sketch<Kmer32bit, u32, FxHasher32> as SeqSketcherT<Kmer32bit> >::Sig >(database_dirpath)?;
                    let bh32 = std::hash::BuildHasherDefault::<FxHasher32>::default();
                    let sketcher = SuperHash2Sketch::<Kmer32bit, u32, FxHasher32>::new(sketch_params, bh32);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, SuperHash2Sketch::<Kmer32bit, u32, FxHasher32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHash2Sketch<Kmer64bit, u64, FxHasher64> as SeqSketcherT<Kmer64bit> >::Sig >(database_dirpath)?;
                    let bh64 = std::hash::BuildHasherDefault::<fxhash::FxHasher64>::default();
                    let sketcher = SuperHash2Sketch::<Kmer64bit, u64, FxHasher64>::new(sketch_params, bh64);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, SuperHash2Sketch::<Kmer64bit, u64, FxHasher64> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHash2Sketch<Kmer16b32bit, u32, FxHasher32> as SeqSketcherT<Kmer16b32bit> >::Sig >(database_dirpath)?;
                    let bh32 = std::hash::BuildHasherDefault::<fxhash::FxHasher32>::default();
                    let sketcher = SuperHash2Sketch::<Kmer16b32bit, u32, FxHasher32>::new(sketch_params, bh32);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, SuperHash2Sketch::<Kmer16b32bit, u32, FxHasher32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search); 
                } 
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));                   
                }          
            }
        }
    } // end match on algo
    //
    return Ok(matcher);  
}  // end of get_sequence_matcher

