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
use crate::dna::dnafiles::{process_file_concat_split, process_file_in_one_block};
use crate::{matcher::*, answer::ReqAnswer};




fn sketch_and_request_dir_compressedkmer<Kmer:CompressedKmerT>(request_dirpath : &Path, filter_params: &FilterParams, 
                    seqdict : &SeqDict, processing_parameters : &ProcessingParams, 
                    hnsw : &Hnsw<Kmer::Val,DistHamming>, knbn : usize, ef_search : usize) -> Matcher
            where Kmer::Val : num::PrimInt + Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
                  KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>, 
                  DistHamming : Distance<Kmer::Val> {
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
    let sketcher = SeqSketcher::new(sketcher_params.get_kmer_size(), sketcher_params.get_sketch_size());
    // create something for likelyhood computation
    let mut matcher = Matcher::new(processing_parameters.get_kmer_size(), sketcher_params.get_sketch_size(), seqdict);
    //
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(1_000);
    // launch process_dir in a thread or async
    crossbeam_utils::thread::scope(|scope| {

        // sequence sending, productor thread
        let mut nb_sent = 0;
        let sender_handle = scope.spawn(move |_|   {
            let res_nb_sent;
            if block_processing {
                res_nb_sent = process_dir(&mut state, &DataType::DNA, request_dirpath, &filter_params, &process_file_in_one_block, &send);
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
                            let signatures = sketcher.sketch_probminhash3a_compressedkmer(&sequencegroup_ref, kmer_hash_fn);                            
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
                            let signatures = sketcher.sketch_probminhash3a_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
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
                     filter_params : &FilterParams, ann_params: &AnnParameters, seqdict : &SeqDict, 
                     nbng : u16, ef_search : usize) -> Result<Matcher, String> {
    //
    let sk_params = processing_params.get_sketching_params();
    log::info!("sketch params reloaded kmer size : {}, sketch size {}", sk_params.get_kmer_size(), sk_params.get_sketch_size());
    //
    let matcher : Matcher;
    // reload hnsw
    log::info!("\n reloading hnsw from {}", database_dirpath.to_str().unwrap());
    if sk_params.get_kmer_size() <= 14 {
        let hnsw = reloadhnsw::reload_hnsw(database_dirpath, ann_params);
        let hnsw = match hnsw {
            Some(hnsw) => hnsw,
            _ => {
                log::error!("\n dna get_sequence_matcher failed to reload hnsw. do you run on DNA data ?");
                panic!("hnsw reload from dump dir {} failed, do you run on DNA data ?", database_dirpath.to_str().unwrap());
            }
        };
        matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit>(&request_dirpath, &filter_params, &seqdict, &processing_params, 
                    &hnsw, nbng as usize, ef_search);
        return Ok(matcher)  
    }
    else if sk_params.get_kmer_size() > 16 {
        let hnsw = reloadhnsw::reload_hnsw(database_dirpath, ann_params);
        let hnsw = match hnsw {
            Some(hnsw) => hnsw,
            _ => {
                log::error!("\n dna get_sequence_matcher failed to reload hnsw. do you run on DNA data ?");
                panic!("hnsw reload from dump dir {} failed, do you run on DNA data ?", database_dirpath.to_str().unwrap());
            }
        };
        matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit>(&request_dirpath, &filter_params, &seqdict, &processing_params,
                    &hnsw, nbng as usize, ef_search);
        return Ok(matcher)  
    }
    else if sk_params.get_kmer_size() == 16 {
        let hnsw = reloadhnsw::reload_hnsw(database_dirpath, ann_params);
        let hnsw = match hnsw {
            Some(hnsw) => hnsw,
            _ => {
                log::error!("\n dna get_sequence_matcher failed to reload hnsw. do you run on DNA data ?");
                panic!("hnsw reload from dump dir {} failed, do you run on DNA data ?", database_dirpath.to_str().unwrap());
            }
        };
        matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit>(&request_dirpath, &filter_params, &seqdict, &processing_params, 
                    &hnsw, nbng as usize, ef_search);
        return Ok(matcher)  
    }
    else {
        log::error!("bad value for kmer size. 15 is not allowed");
        Err(String::from("bad value for kmer size"))
    }
}  // end of get_sequence_matcher

