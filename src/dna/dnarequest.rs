//! gathers dna request stuf
//! 
//! 

// TODO: optimize nb_bases_threshold



use std::path::{Path, PathBuf};
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use crossbeam::sync::Parker;

use std::time::Duration;
use std::time::SystemTime;
use cpu_time::ProcessTime;

// for parallellism
use crossbeam::queue::ArrayQueue;

use serde::{Serialize, de::DeserializeOwned};
use std::fmt::Debug;

#[allow(unused)]
use fxhash::*;

use hnsw_rs::prelude::*;

use kmerutils::base::{kmergenerator::*, Kmer32bit, Kmer64bit, CompressedKmerT};
use kmerutils::sketching::setsketchert::*;

use probminhash::setsketcher::SetSketchParams;

use crate::utils::*;
use crate::dna::dnafiles::{process_file_by_sequence, process_buffer_by_sequence, process_buffer_in_one_block, process_file_in_one_block};
use crate::{matcher::*, answer::ReqAnswer};




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
                computing_params : &ComputingParams, 
                hnsw : &Hnsw< <Sketcher as SeqSketcherT<Kmer>>::Sig ,DistHamming>, 
                knbn : usize, ef_search : usize) -> Matcher

            where   Kmer::Val : num::PrimInt + Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
                    Sketcher : SeqSketcherT<Kmer> + Clone + Send + Sync,
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
    let outpath = PathBuf::from(outname);
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
    // but if request files are very large we must not have too many files in memory
    // pario is used also to limit the number of file loaded simultanuously.
    let request_block_size = match computing_params.get_parallel_io() {
        true => { 2000.min(computing_params.get_nb_files_par()) },
        _    => { 20 },
    };
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let pool_nb_thread = pool.current_num_threads();
    let mut nb_max_threads = computing_params.get_sketching_nbthread();
    if nb_max_threads == 0 {
        nb_max_threads = request_block_size.min(pool_nb_thread).max(1);
    }
    log::info!("nb threads in pool : {:?}, using nb threads : {}", pool_nb_thread, nb_max_threads);
    log::debug!("request_block_size : {}", request_block_size);
    //
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
    let mut nb_received = 0;
    //
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(request_block_size);
    // to send sketch result to a collector task
    let (request_sender , request_receiver) = 
        crossbeam_channel::bounded::<RequestMsg<<Sketcher as SeqSketcherT<Kmer>>::Sig>>(request_block_size+1);
    //
    // launch process_dir in a thread or async
    //
    pool.scope(|scope| {
        // sequence sending, productor thread
        scope.spawn(|_|   {
            let res_nb_sent;
            match block_processing {
                true => {
                    log::info!("dnarequest::sketchandstore_dir_compressedkmer : block processing");
                    if computing_params.get_parallel_io() {
                        let nb_files_by_group = computing_params.get_nb_files_par();
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
                    if computing_params.get_parallel_io() {
                        let nb_files_by_group = computing_params.get_nb_files_par();
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

        // sequence reception, sketching/request consumer thread. The real hnsw request is delegated to a special thread
        scope.spawn( |scope| {
            // we must read messages, sketch and insert into hnsw
            let parker: Parker = Parker::new();
            // we can create a new thread for at least nb_bases_thread_threshold bases.
            let nb_bases_thread_threshold : usize = 100_000_000;
            // a bounded blocking channel to limit the number of threads to pool_nb_thread.
            // at thread creation we send a msg into queue, at thread end we receive a msg.
            // So the length of the channel is number of active thread
            let (thread_token_sender, thread_token_receiver) = crossbeam_channel::bounded::<u32>(nb_max_threads);
            //
            let sketcher_queue = ArrayQueue::<Vec<IdSeq>>::new(request_block_size);
            let mut read_more = true;
            let mut nb_base_in_queue = 0;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                if !sketcher_queue.is_full() {
                    let res_receive = receive.recv();
                    match res_receive {
                        Err(_) => { read_more = false;
                            log::debug!("end of sequence request reception");
                        },
                        Ok(idsequences) => {
                            // concat the new idsketch in sketcher queue.
                            log::debug!("request reception nb seq : {}", idsequences.len());
                            nb_base_in_queue += idsequences.iter().fold(0, |acc, s| acc + s.get_seq_len());
                            log::debug!("nb_base_in_queue : {}", nb_base_in_queue);
                            let p_res = sketcher_queue.push(idsequences);
                            if p_res.is_err() {
                                log::error!("aborting in thread dispatcher, cannot push in sketcher queue");
                                std::process::abort();
                            };
                        }
                    };
                }
                //
                if read_more == false && sketcher_queue.len() == 0 && thread_token_sender.len() == 0 {
                    log::info!("no more read to come, no more data to sketch, no more sketcher running ...");
                    break;
                }
                else if read_more == false && sketcher_queue.len() == 0 {
                    // TODO here we just have to wait all thread end, should use a condvar, now we check every second
                    parker.park_timeout(Duration::from_millis(1000));
                }
                //
                // if request_queue is beyond threshold size we can go to threaded sketching and threading insertion
                // we will spawn a thread that must do sketching of request and send request singature to hnsw dedicated thread
                //
                if nb_base_in_queue >= nb_bases_thread_threshold || sketcher_queue.is_full() || (!read_more && !sketcher_queue.is_empty()) {
                    let nb_free_thread = thread_token_sender.capacity().unwrap() as i32 - thread_token_sender.len() as i32;
                    log::debug!("nb_free_thread : {}",nb_free_thread);
                    let request_sender_clone = request_sender.clone();
                    let sketcher_clone = sketcher.clone();
                    // transfer from  sketcher_queue to a local queue that will be move to spawn task
                    let q_len = sketcher_queue.len();
                    assert!(q_len > 0);
                    let to_pop = q_len;
                    log::debug!("popping to local_queue {}, sketching_queue.len() : {}", to_pop, sketcher_queue.len());
                    let mut local_queue = Vec::<Vec<IdSeq> >::with_capacity(to_pop);
                    for _ in 0..to_pop {
                        let res_pop = sketcher_queue.pop();
                        let mut nb_popped_bases : usize = 0;
                        if res_pop.is_some() {
                            nb_popped_bases += res_pop.as_ref().unwrap().iter().fold(0, |acc, s| acc + s.get_seq_len());
                            local_queue.push(res_pop.unwrap());
                        }
                        assert!(nb_base_in_queue >= nb_popped_bases);
                        nb_base_in_queue -= nb_popped_bases; 
                    }
                    log::debug!("after popping , nb_base_in_queue : {}", nb_base_in_queue);
                    let q_len = local_queue.len();
                    assert!(q_len > 0);
                    //
                    let _ = thread_token_sender.send(1);
                    log::debug!("sending token sender tokens = {:?}, nb_token_reciever : {}", thread_token_sender.len(), thread_token_receiver.len());
                    let thread_token_receiver_cloned = thread_token_receiver.clone();
                    //
                    scope.spawn(move |_|  {
                        log::trace!("spawning thread on nb files : {}", local_queue.len());
                        match block_processing { 
                            true => {
                                for v in &local_queue {
                                    assert_eq!(v.len(), 1);
                                }
                                let sequencegroup_ref : Vec<&Sequence> = local_queue.iter().map(|s| s[0].get_sequence_dna().unwrap()).collect();
                                // collect Id
                                let seq_item : Vec<ItemDict> = local_queue.iter().map(|s| ItemDict::new(Id::new(s[0].get_path(), s[0].get_fasta_id()), s[0].get_seq_len())).collect();
                                // computes hash signature
                                let signatures = sketcher_clone.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
                                // now we must send signatures, rak and seq_items to request collector
                                for i in 0..signatures.len() {
                                    let request_msg = RequestMsg::new(signatures[i].clone() , seq_item[i].clone());
                                    let _ = request_sender_clone.send(request_msg);
                                }
                            },
                            //
                            false => { // means we are in seq by seq mode inside a file. We get one signature by msg (or file)
                                for i in 0..local_queue.len() {
                                    let sequencegroup_ref : Vec<&Sequence> = local_queue[i].iter().map(|s| s.get_sequence_dna().unwrap()).collect();
                                    let seq_len = sequencegroup_ref.iter().fold(0, | acc, s | acc + s.size());
                                    // collect Id
                                    let seq_item : ItemDict = ItemDict::new(Id::new(local_queue[i][0].get_path(), local_queue[i][0].get_fasta_id()), seq_len);
                                    // computes hash signature
                                    let signatures = sketcher_clone.sketch_compressedkmer_seqs(&sequencegroup_ref, kmer_hash_fn);
                                    log::debug!("msg num : {}, nb seq : {}, path : {:?}", i, local_queue[i].len(), seq_item.get_id().get_path());
                                    // now we must send signatures, rank and seq_items to request collector
                                    assert_eq!(signatures.len(), 1);
                                    let request_msg = RequestMsg::new(signatures[0].clone() , seq_item.clone());
                                    let _ = request_sender_clone.send(request_msg);
                                }
                            },  // end by seq case
                        }; // end match
                        // cleaning
                        while let Some(v) = local_queue.pop() {
                            drop(v);
                        }
                        local_queue.clear();  
                        drop(local_queue); 
                        // we free a token in queue
                        let res = thread_token_receiver_cloned.recv();
                        log::debug!("thread_token_receiver_cloned.len = {:?}", thread_token_receiver_cloned.len());
                        match res {
                            Ok(_)   => { log::debug!("thread ending OK");},
                            Err(_)  => { log::error!("thread exit with a token error, aborting");
                                        std::process::abort();
                                    },
                        }                         
                    }); // end spawned thread
                }
            } // end while 
            //
            drop(receive);
            drop(request_sender);
        }); // end of receptor sketcher thread
        //
        // now we spawn a new thread that is devoted to is devoted to hnsw search
        //
        scope.spawn(|_| {
            // Sig is basic item of a signature , VecSig is a vector of such items
            type VecSig<Sketcher, Kmer>  = Vec< <Sketcher as SeqSketcherT<Kmer>>::Sig>;
            let nb_request_guess = request_block_size.min(1000);
            let mut request_store = Vec::<VecSig<Sketcher,Kmer>>::with_capacity(nb_request_guess);
            let mut itemv = Vec::<ItemDict>::with_capacity(nb_request_guess);
            //
            let mut read_more = true;
            let mut nb_request = 0;
            let mut nb_answer = 0;

            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = request_receiver.recv();
                match res_receive {
                    Err(_recv_error) =>  { read_more = false;
                                                    log::debug!("end of request reception");
                    }
                    Ok(to_insert) => {
                        request_store.push(to_insert.sketch_and_rank);
                        itemv.push(to_insert.item);
                        nb_request += 1;
                        log::debug!("request collector received nb_received : {}", nb_received);
                    }
                }
                // we do all request at end beccause of scheduling pb between par_iter and explicit thread in rayon 
                if read_more == false {
                    log::debug!("recieved requests nb : {:?}", request_store.len());
                    let mut data_for_hnsw = Vec::<VecSig<Sketcher,Kmer> >::with_capacity(request_store.len());
                    for i in 0..request_store.len() {
                        data_for_hnsw.push(request_store[i].clone());
                    }
                    let knn_neighbours = hnsw.parallel_search(&data_for_hnsw, knbn, ef_search);  
                    // construct and dump answers
                    for i in 0..knn_neighbours.len() {
                        let answer = ReqAnswer::new(nb_answer+i, itemv[i].clone(), &knn_neighbours[i]);
                        if answer.dump(&seqdict, out_threshold, &mut outfile).is_err() {
                            log::info!("could not dump answer for request id {}", answer.get_request_id().get_id().get_fasta_id());
                        }
                        // store in matcher. remind that each i corresponds to a request
                        let candidates = knn_neighbours[i].iter().map(|n| SequenceMatch::new(seqdict.0[n.d_id].clone(), n.get_distance())).collect();
                        matcher.insert_sequence_match(itemv[i].clone(), candidates);
                    }
                    nb_answer += knn_neighbours.len();
                    let _ = outfile.flush();
                    request_store.clear();
                    itemv.clear();
                }
            }
            let _ = outfile.flush();            
            // transfer nb_received to global nb_received
            nb_received = nb_request;
            log::info!("request collector , exiting after answering to nb_request {}", nb_received);
        }); // end of collector thread
    });  // end of pool.scope
    //
    log::debug!("sketch_and_request_dir_compressedkmer, nb_sent = {}, nb_received = {}", nb_sent, nb_received);
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
    //
    let hnswio_res = reloadhnsw::get_hnswio(database_dirpath);
    if hnswio_res.is_err() {
        std::panic!("error : {:?}", hnswio_res.err());
    }
    let mut hnswio = hnswio_res.unwrap();
    // allocate sketcher and get matcher
    match sketch_params.get_algo() {
        SketchAlgo::PROB3A => {
            match sketch_params.get_kmer_size() {
                1..=14 => {
                    let hnsw_res = hnswio.load_hnsw::< <Kmer32bit as CompressedKmerT>::Val, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();
                    let sketcher = ProbHash3aSketch::<Kmer32bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, ProbHash3aSketch::<Kmer32bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &computing_params, 
                            &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw_res = hnswio.load_hnsw::< <Kmer64bit as CompressedKmerT>::Val, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();
                    let sketcher = ProbHash3aSketch::<Kmer64bit>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, ProbHash3aSketch::<Kmer64bit> >(&request_dirpath, sketcher, 
                            &filter_params, &seqdict, &processing_params, &computing_params,
                            &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw_res = hnswio.load_hnsw::< <Kmer16b32bit as CompressedKmerT>::Val, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();
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
                    let hnsw_res = hnswio.load_hnsw::< <SuperHashSketch<Kmer32bit, f32> as SeqSketcherT<Kmer32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();
                    let sketcher = SuperHashSketch::<Kmer32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, SuperHashSketch::<Kmer32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw_res = hnswio.load_hnsw::< <SuperHashSketch<Kmer64bit, f32> as SeqSketcherT<Kmer64bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();                    
                    let sketcher = SuperHashSketch::<Kmer64bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, SuperHashSketch::<Kmer64bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw_res = hnswio.load_hnsw::< <SuperHashSketch<Kmer16b32bit, f32> as SeqSketcherT<Kmer16b32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();  
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
        //
        SketchAlgo::OPTDENS => {
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw_res = hnswio.load_hnsw::< <OptDensHashSketch<Kmer32bit, f32> as SeqSketcherT<Kmer32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();
                    let sketcher = OptDensHashSketch::<Kmer32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, OptDensHashSketch::<Kmer32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw_res = hnswio.load_hnsw::< <OptDensHashSketch<Kmer64bit, f32> as SeqSketcherT<Kmer64bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();                    
                    let sketcher = OptDensHashSketch::<Kmer64bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, OptDensHashSketch::<Kmer64bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw_res = hnswio.load_hnsw::< <OptDensHashSketch<Kmer16b32bit, f32> as SeqSketcherT<Kmer16b32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();  
                    let sketcher = OptDensHashSketch::<Kmer16b32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, OptDensHashSketch::<Kmer16b32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search); 
                } 
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));                   
                }  
            }
        } // end case OptDensHashSketch
        //
        SketchAlgo::HLL => {
            // 16 bits is sufficient up to kmer of size 30!
            let mut hll_params = SetSketchParams::default();
            if hll_params.get_m() < sketch_params.get_sketch_size() as u64 {
                log::warn!("!!!!!!!!!!!!!!!!!!!  need to adjust hll parameters!");
            }
            hll_params.set_m(sketch_params.get_sketch_size());
            //
            let nb_cpus = num_cpus::get();
            log::info!("nb cpus : {}", nb_cpus);
            let nb_iter_thtreads = ((nb_cpus.max(4) - 4) / computing_params.get_sketching_nbthread().max(1)).max(1);
            let hll_seqs_threading = HllSeqsThreading::new(nb_iter_thtreads, 10_000_000);
            log::info!("HllSeqsThreading : {:?}", hll_seqs_threading);
            //
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw_res = hnswio.load_hnsw::< <HyperLogLogSketch<Kmer32bit, u16> as SeqSketcherT<Kmer32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();  
                    let sketcher = HyperLogLogSketch::<Kmer32bit, u16>::new(sketch_params, hll_params, hll_seqs_threading);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, HyperLogLogSketch::<Kmer32bit, u16> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw_res = hnswio.load_hnsw::< <HyperLogLogSketch<Kmer64bit, u16> as SeqSketcherT<Kmer64bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();  
                    let sketcher = HyperLogLogSketch::<Kmer64bit, u16>::new(sketch_params, hll_params, hll_seqs_threading);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, HyperLogLogSketch::<Kmer64bit, u16> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw_res = hnswio.load_hnsw::< <HyperLogLogSketch<Kmer16b32bit, u16> as SeqSketcherT<Kmer16b32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap(); 
                    let sketcher = HyperLogLogSketch::<Kmer16b32bit, u16>::new(sketch_params, hll_params, hll_seqs_threading);
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
                    let hnsw_res = hnswio.load_hnsw::< <SuperHash2Sketch<Kmer32bit, u32, FxHasher32> as SeqSketcherT<Kmer32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap(); 
                    let bh32 = std::hash::BuildHasherDefault::<fxhash::FxHasher32>::default();
                    let sketcher = SuperHash2Sketch::<Kmer32bit, u32, FxHasher32>::new(sketch_params, bh32);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, SuperHash2Sketch::<Kmer32bit, u32, FxHasher32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw_res = hnswio.load_hnsw::< <SuperHash2Sketch<Kmer64bit, u64, FxHasher64> as SeqSketcherT<Kmer64bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap(); 
                    let bh64 = std::hash::BuildHasherDefault::<fxhash::FxHasher64>::default();
                    let sketcher = SuperHash2Sketch::<Kmer64bit, u64, FxHasher64>::new(sketch_params, bh64);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, SuperHash2Sketch::<Kmer64bit, u64, FxHasher64> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw_res = hnswio.load_hnsw::< <SuperHash2Sketch<Kmer16b32bit, u32, FxHasher32> as SeqSketcherT<Kmer16b32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap(); 
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
        } // end superminhash2
        //
        SketchAlgo::REVOPTDENS => {
            match sketch_params.get_kmer_size() {
                0..=14 => {
                    let hnsw_res = hnswio.load_hnsw::< <RevOptDensHashSketch<Kmer32bit, f32> as SeqSketcherT<Kmer32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();
                    let sketcher = RevOptDensHashSketch::<Kmer32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer32bit, RevOptDensHashSketch::<Kmer32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                17..=32 => {
                    let hnsw_res = hnswio.load_hnsw::< <RevOptDensHashSketch<Kmer64bit, f32> as SeqSketcherT<Kmer64bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();                    
                    let sketcher = RevOptDensHashSketch::<Kmer64bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer64bit, RevOptDensHashSketch::<Kmer64bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search);
                }
                16 => {
                    let hnsw_res = hnswio.load_hnsw::< <RevOptDensHashSketch<Kmer16b32bit, f32> as SeqSketcherT<Kmer16b32bit> >::Sig, DistHamming>();
                    if hnsw_res.is_err() {
                        std::panic!("error : {:?}", hnsw_res.err());
                    };
                    let hnsw = hnsw_res.unwrap();  
                    let sketcher = RevOptDensHashSketch::<Kmer16b32bit, f32>::new(sketch_params);
                    matcher = sketch_and_request_dir_compressedkmer::<Kmer16b32bit, RevOptDensHashSketch::<Kmer16b32bit, f32> >(&request_dirpath, sketcher, 
                        &filter_params, &seqdict, &processing_params, &computing_params,
                        &hnsw, nbng as usize, ef_search); 
                } 
                _ => {
                    log::error!("bad value for kmer size. 15 is not allowed");
                    return Err(String::from("bad value for kmer size"));                   
                }  
            }
        } // end case RevOptDensHashSketch
        
    } // end match on algo
    //
    return Ok(matcher);  
}  // end of get_sequence_matcher

