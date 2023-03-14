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

use hnsw_rs::prelude::*;

use kmerutils::base::{kmergenerator::*, Kmer32bit, Kmer64bit, CompressedKmerT};
use kmerutils::sketching::seqsketchjaccard::*;

use crate::utils::*;
use crate::dna::dnafiles::{process_file_concat_split, process_buffer_in_one_block, process_file_in_one_block};
use crate::{matcher::*, answer::ReqAnswer};




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
    let outname = "gsearch.answers";
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
    let mut request_queue : Vec<IdSeq>= Vec::with_capacity(request_block_size);
    let mut state = ProcessingState::new();
    
    //
    // Sketcher allocation, we do not need reverse complement hashing as we sketch assembled genomes. (Jianshu Zhao)
    //
    let kmer_hash_fn = | kmer : &Kmer | -> Kmer::Val {
        let mask : Kmer::Val = num::NumCast::from::<u64>((0b1 << 2*kmer.get_nb_base()) - 1).unwrap();
        let hashval = kmer.get_compressed_value() & mask;
        hashval
    };
    // create something for likelyhood computation
    let mut matcher = Matcher::new(processing_parameters.get_kmer_size(), sketcher_params.get_sketch_size(), seqdict);
    //
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(request_block_size + 10);
    // launch process_dir in a thread or async
    crossbeam_utils::thread::scope(|scope| {

        // sequence sending, productor thread
        let mut nb_sent = 0;
        let sender_handle = scope.spawn(move |_|   {
            let res_nb_sent;
            if block_processing {
                if other_params.get_parallel_io() {
                    let nb_files_by_group = other_params.get_nb_files_par();
                    log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir_parallel, nb_files in parallel : {}", nb_files_by_group);
                    res_nb_sent = process_dir_parallel(&mut state, &DataType::DNA, request_dirpath, filter_params, 
                                    nb_files_by_group, &process_buffer_in_one_block, &send);
                } // end case parallel io
                else {
                    res_nb_sent = process_dir(&mut state, &DataType::DNA, request_dirpath, &filter_params, &process_file_in_one_block, &send);
                }
            }
            else {
                log::info!("processing by concat and split");
                res_nb_sent = process_dir(&mut state, &DataType::DNA,  request_dirpath, &filter_params, &process_file_concat_split, &send);
            }
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
            Box::new(nb_sent)
        });
        // sequence reception, consumer thread
        let receptor_handle = scope.spawn( |_| {
            let mut nb_request = 0;
            // we must read messages, sketch and insert into hnsw
            let mut read_more = true;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = receive.recv();
                match res_receive {
                    Err(_) => { read_more = false;
                        // sketch the content of  insertion_queue if not empty
                        if request_queue.len() > 0 {
                            let sequencegroup_ref : Vec<&Sequence> = request_queue.iter().map(|s| s.get_sequence_dna().unwrap()).collect();                   
                            let seq_item : Vec<ItemDict> = request_queue.iter().map(|s| ItemDict::new(Id::new(s.get_path(), s.get_fasta_id()), s.get_seq_len())).collect();
                            if log::log_enabled!(log::Level::Debug) {
                                for s in &seq_item  {
                                    log::debug!("treating request : {}", s.get_id().get_path());
                                }
                            }   
                            let signatures = sketcher.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);                            
                            // we have Vec<u64> signatures we must go back to a vector of IdSketch for hnsw insertion
                            let knn_neighbours  = hnsw.parallel_search(&signatures, knbn, ef_search);
                            for i in 0..knn_neighbours.len() {
                                let answer = ReqAnswer::new(nb_request+i, seq_item[i].clone(), &knn_neighbours[i]);
                                if answer.dump(&seqdict, out_threshold, &mut outfile).is_err() {
                                    log::info!("could not dump answer for request id {}", answer.get_request_id().get_id().get_fasta_id());
                                }
                                // store in matcher. remind that each i corresponds to a request
                                let candidates = knn_neighbours[i].iter().map(|n| SequenceMatch::new(seqdict.0[n.d_id].clone(), n.get_distance())).collect();
                                matcher.insert_sequence_match(seq_item[i].clone(), candidates);
                            }
                            //  dump results
                            nb_request += signatures.len();
                        }
                    }
                    Ok(mut idsequences) => {
                        // concat the new idsketch in insertion queue.
                        request_queue.append(&mut idsequences);
                        // if request_queue is beyond threshold size we can go to threaded sketching and threading insertion
                        if request_queue.len() > request_block_size {
                            let sequencegroup_ref : Vec<&Sequence> = request_queue.iter().map(|s| s.get_sequence_dna().unwrap()).collect();
                            // collect Id
                            let seq_item : Vec<ItemDict> = request_queue.iter().map(|s| ItemDict::new(Id::new(s.get_path(), s.get_fasta_id()), s.get_seq_len())).collect();
                            // computes hash signature
                            let signatures = sketcher.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
                            // we have Vec<u64> signatures we must go back to a vector of IdSketch, inserting unique id, for hnsw insertion
                            // parallel search
                            let knn_neighbours  = hnsw.parallel_search(&signatures, knbn, ef_search);
                            // construct and dump answers
                            for i in 0..knn_neighbours.len() {
                                let answer = ReqAnswer::new(nb_request+i, seq_item[i].clone(), &knn_neighbours[i]);
                                if answer.dump(&seqdict, out_threshold, &mut outfile).is_err() {
                                    log::info!("could not dump answer for request id {}", answer.get_request_id().get_id().get_fasta_id());
                                }
                                // store in matcher. remind that each i corresponds to a request
                                let candidates = knn_neighbours[i].iter().map(|n| SequenceMatch::new(seqdict.0[n.d_id].clone(), n.get_distance())).collect();
                                matcher.insert_sequence_match(seq_item[i].clone(), candidates);
                            }
                            nb_request += signatures.len();
                            request_queue.clear();
                        }
                    }
                }
            } // end while 
            //
            Box::new(nb_request)
        }); // end of receptor thread
        // now we must join handles
        let nb_sent = sender_handle.join().unwrap();
        let nb_received = receptor_handle.join().unwrap();
        log::debug!("sketch_and_request_dir, nb_sent = {}, nb_received = {}", nb_sent, nb_received);
        if nb_sent != nb_received {
            log::error!("an error occurred  nb msg sent : {}, nb msg received : {}", nb_sent, nb_received);
        }
    }  // end of closure in scope
    ).unwrap();  // end of scope
    //
    log::info!("matcher collected {} answers", matcher.get_nb_sequence_match());
    let cpu_time = cpu_start.elapsed().as_secs();
    log::info!("process_dir : cpu time(s) {}", cpu_time);
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
    log::info!("process_dir : elapsed time(s) {}", elapsed_t);
    //
    matcher
} // end of sketch_and_request_dir_compressedkmer 





// This function returns paired sequence by probminhash and hnsw 
pub fn get_sequence_matcher(request_dirpath : &Path, database_dirpath : &Path, processing_params : &ProcessingParams,
                     filter_params : &FilterParams, ann_params: &AnnParameters, other_params : &ComputingParams, seqdict : &SeqDict, 
                     nbng : u16, ef_search : usize) -> Result<Matcher, String> {
    //
    let sketch_params = processing_params.get_sketching_params();
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
                    let hnsw = reloadhnsw::reload_hnsw::< <Kmer32bit as CompressedKmerT>::Val>(database_dirpath, ann_params)?;
                    let sketcher = ProbHash3aSketch::<Kmer32bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, ProbHash3aSketch::<Kmer32bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &other_params, 
                            &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <Kmer64bit as CompressedKmerT>::Val>(database_dirpath, ann_params)?;
                    let sketcher = ProbHash3aSketch::<Kmer64bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, ProbHash3aSketch::<Kmer64bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &other_params,
                            &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <Kmer16b32bit as CompressedKmerT>::Val>(database_dirpath, ann_params)?;
                    let sketcher = ProbHash3aSketch::<Kmer16b32bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, ProbHash3aSketch::<Kmer16b32bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &other_params,
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
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHashSketch<Kmer32bit, f32> as SeqSketcherT<Kmer32bit> >::Sig >(database_dirpath, ann_params)?;
                    let sketcher = SuperHashSketch::<Kmer32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, SuperHashSketch::<Kmer32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &other_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHashSketch<Kmer64bit, f32> as SeqSketcherT<Kmer64bit> >::Sig >(database_dirpath, ann_params)?;
                    let sketcher = SuperHashSketch::<Kmer64bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, SuperHashSketch::<Kmer64bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &other_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHashSketch<Kmer16b32bit, f32> as SeqSketcherT<Kmer16b32bit> >::Sig >(database_dirpath, ann_params)?;
                    let sketcher = SuperHashSketch::<Kmer16b32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, SuperHashSketch::<Kmer16b32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &other_params,
                        &hnsw, nbng as usize, ef_search); 
                } 
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));                   
                }          
            }
        }  // end match SUPER
        SketchAlgo::SUPER2 => {
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHash2Sketch<Kmer32bit, u32> as SeqSketcherT<Kmer32bit> >::Sig >(database_dirpath, ann_params)?;
                    let sketcher = SuperHash2Sketch::<Kmer32bit, u32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, SuperHash2Sketch::<Kmer32bit, u32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &other_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHash2Sketch<Kmer64bit, u64> as SeqSketcherT<Kmer64bit> >::Sig >(database_dirpath, ann_params)?;
                    let sketcher = SuperHash2Sketch::<Kmer64bit, u64>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, SuperHash2Sketch::<Kmer64bit, u64> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &other_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw = reloadhnsw::reload_hnsw::< <SuperHash2Sketch<Kmer16b32bit, u32> as SeqSketcherT<Kmer16b32bit> >::Sig >(database_dirpath, ann_params)?;
                    let sketcher = SuperHash2Sketch::<Kmer16b32bit, u32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, SuperHash2Sketch::<Kmer16b32bit, u32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &other_params,
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

