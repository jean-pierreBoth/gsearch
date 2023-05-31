//! amino acid (AA) sketching


use std::time::{SystemTime};
use cpu_time::ProcessTime;

// for multithreading
use std::sync::Arc;
use crossbeam_channel::*;
use concurrent_queue::{ConcurrentQueue, PushError};

use serde::{de::DeserializeOwned, Serialize};

use std::fmt::{Debug};
use std::path::{PathBuf};
use cpu_time::ThreadTime;


// our crate
use hnsw_rs::prelude::*;

use probminhash::{setsketcher::SetSketchParams};

use kmerutils::base::kmertraits::*;
use kmerutils::aautils::kmeraa::*;
use kmerutils::aautils::seqsketchjaccard::{SeqSketcherAAT, ProbHash3aSketch, SuperHashSketch, HyperLogLogSketch};
use kmerutils::sketcharg::DataType;

use crate::utils::{idsketch::*, reloadhnsw};
use crate::utils::dumpload::*;

use crate::utils::files::{process_dir, process_dir_parallel, ProcessingState};
use crate::aa::aafiles::{process_aafile_in_one_block, process_aabuffer_in_one_block, process_aabuffer_by_sequence, process_aafile_by_sequence};


use crate::utils::parameters::*;

// Sig is basic item of a signature , VecSig is a vector of such items
type VecSig<Sketcher, Kmer>  = Vec< <Sketcher as SeqSketcherAAT<Kmer>>::Sig>;

// a type to describe msessage to collector task
struct CollectMsg<Sig> {
    // signature and sequence rank
    pub skecth_and_rank : (Vec<Sig>, usize),
    // id of sequence
    pub item : ItemDict,
}


impl <Sig> CollectMsg<Sig> {

    pub fn new(skecth_and_rank : (Vec<Sig>, usize), item : ItemDict) -> Self {
        CollectMsg{skecth_and_rank, item}
    }

} // end of CollectMsg



fn sketchandstore_dir_compressedkmer<Kmer:CompressedKmerT+KmerBuilder<Kmer>, Sketcher : SeqSketcherAAT<Kmer>>(hnsw_pb : &PathBuf, sketcher : Sketcher ,
                    filter_params: &FilterParams, 
                    processing_params : &ProcessingParams, 
                    other_params : &ComputingParams) 

        where   Sketcher : SeqSketcherAAT<Kmer> + Clone + Send + Sync,
                <Sketcher as SeqSketcherAAT<Kmer>>::Sig : 'static + Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>, 
                DistHamming : Distance< <Sketcher as  SeqSketcherAAT<Kmer>>::Sig>,
                {
    //
    //
    log::info!("sketchandstore_dir {}", hnsw_pb.to_str().unwrap());
    let start_t = SystemTime::now();
    let cpu_start: ProcessTime = ProcessTime::now();
    //
    let block_processing = processing_params.get_block_flag();
    // a queue of signature waiting to be inserted , size must be sufficient to benefit from threaded probminhash and insert
    // and not too large to spare memory. If parallel_io is set dimension message queue to size of group
    // for files of size more than Gb we must use pario to limit memory, but leave enough msg in queue to get // sketch and insertion 
    let insertion_block_size = match other_params.get_parallel_io() {
        true => { 2000.min(1 + other_params.get_nb_files_par()) },
        _    => { 2000 },
    };
    log::debug!("insertion_block_size : {}", insertion_block_size);
    //
    let mut hnsw : Hnsw::< <Sketcher as SeqSketcherAAT<Kmer>>::Sig, DistHamming>;
    let mut state :ProcessingState;
    //
    let toprocess_path : PathBuf;
    //
    if other_params.get_adding_mode() {
        let toprocess_str = other_params.get_add_dir();
        toprocess_path = PathBuf::from(toprocess_str);
        // in this case we must reload
        log::info!("dnasketch::sketchandstore_dir_compressedkmer will add new data from directory {:?} to hnsw dir {:?}", toprocess_path, hnsw_pb);
        let hnsw_opt = reloadhnsw::reload_hnsw(&hnsw_pb);
        if hnsw_opt.is_err() {
            log::error!("cannot reload hnsw from directory : {:?}", &hnsw_pb);
            std::process::exit(1);
        }
        else { 
            hnsw = hnsw_opt.unwrap();
        }
        let reload_res = ProcessingState::reload_json(&hnsw_pb);
        if reload_res.is_ok() {
            state = reload_res.unwrap();
        }
        else {
            log::error!("dnasketch::cannot reload processing state (file 'processing_state.json' from directory : {:?}", &hnsw_pb);
            std::process::exit(1);           
        }
    } 
    else {
        // creation mode
        toprocess_path = hnsw_pb.clone();
        let hnsw_params = processing_params.get_hnsw_params();
        hnsw = Hnsw::< <Sketcher as SeqSketcherAAT<Kmer>>::Sig, DistHamming>::new(hnsw_params.get_max_nb_connection() as usize , hnsw_params.capacity , 16, hnsw_params.get_ef(), DistHamming{});
        state = ProcessingState::new();
    }
    //
    let mut seqdict : SeqDict;
    seqdict = get_seqdict(&hnsw_pb, other_params).unwrap();
    //
    // where do we dump hnsw* seqdict and so on
    // If in add mode we dump where is already an hnsw database
    // If creation mode we dump in .
    //
    let dump_path= if other_params.get_adding_mode() {
        hnsw_pb.clone()
    } else {
        PathBuf::from(".")
    };
    let dump_path_ref = &dump_path;
    //
    hnsw.set_extend_candidates(true);
    hnsw.set_keeping_pruned(false);
    //
    // Sketcher allocation, we  need reverse complement
    //
    let kmer_hash_fn = | kmer : &Kmer | -> Kmer::Val {
        let canonical =  kmer.reverse_complement().min(*kmer);
        let mask : Kmer::Val = num::NumCast::from::<u64>((0b1 << 2*kmer.get_nb_base()) - 1).unwrap();
        let hashval = canonical.get_compressed_value() & mask;
        hashval
    };
    //
    let mut nb_sent : usize = 0;       // this variable is moved in sender  io thread
    let mut nb_received : usize = 0;   // this variable is moved in collector thread!
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(insertion_block_size+1);
    // to send sketch result to a collector task
    let (collect_sender , collect_receiver) = 
            crossbeam_channel::bounded::<CollectMsg<<Sketcher as SeqSketcherAAT<Kmer>>::Sig>>(insertion_block_size+1);
    //
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let pool_nb_thread = pool.current_num_threads();
    let nb_max_threads = insertion_block_size.min(pool_nb_thread);
    log::info!("nb threads in pool : {:?}, using nb threads : {}", pool_nb_thread, nb_max_threads);

    // launch process_dir in a thread or async
    pool.scope(|scope| {
        // sequence sending, productor thread
        scope.spawn( |_|   {
            let start_t_prod = SystemTime::now();
            let cpu_start_prod = ThreadTime::now();
            let res_nb_sent;
            match block_processing {
                true => {
                    log::info!("dnasketch::sketchandstore_dir_compressedkmer : block processing");
                    if other_params.get_parallel_io() {
                        let nb_files_by_group = other_params.get_nb_files_par();
                        log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir_parallel, nb_files in parallel : {}", nb_files_by_group);
                        res_nb_sent = process_dir_parallel(&mut state, &DataType::DNA,  &toprocess_path, filter_params, 
                                        nb_files_by_group, &process_aabuffer_in_one_block, &send);
                    } // end case parallel io
                    else {
                        log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir serial");
                        res_nb_sent = process_dir(&mut state, &DataType::DNA, &toprocess_path, filter_params, 
                                        &process_aafile_in_one_block, &send);
                    }                
                }
                false => {
                    log::info!("dnasketch::sketchandstore_dir_compressedkmer : seq by seq processing");
                    if other_params.get_parallel_io() {
                        let nb_files_by_group = other_params.get_nb_files_par();
                        log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir_parallel, nb_files in parallel : {}", nb_files_by_group);
                        res_nb_sent = process_dir_parallel(&mut state, &DataType::DNA,  &toprocess_path, filter_params, 
                                nb_files_by_group, &process_aabuffer_by_sequence, &send);
                    }
                    else {
                        log::info!("processing by sequence, sketching whole file globally");
                        res_nb_sent = process_dir(&mut state, &DataType::DNA, &toprocess_path, filter_params,
                                     &process_aafile_by_sequence, &send);
                    }
                }                
            };
            match res_nb_sent {
                Ok(nb_really_sent) => {
                    nb_sent = nb_really_sent;
                    log::info!("process_dir processed nb files (msg) : {}", nb_sent);
                    println!("process_dir processed nb files (msg): {}", nb_sent);
                }
                Err(_) => {
                    nb_sent = 0;
                    println!("some error occured in process_dir");
                    log::error!("some error occured in process_dir, aborting in thread sending files/sequences");
                    std::process::abort();
                }
            };
            drop(send);
            state.elapsed_t =  start_t_prod.elapsed().unwrap().as_secs() as f32;
            log::info!("sender processed in  system time(s) : {}, cpu time(s) : {}", state.elapsed_t, cpu_start_prod.elapsed().as_secs());
            // dump processing state in the current directory
            let _ = state.dump_json(dump_path_ref);
        });
        //
        // sequence reception, consumer thread
        //
        scope.spawn( |scope| {
            // we need a bounded queue (no need for a concurrent queue yet) to be able to block the thread if max number of thread is reached
            let insertion_queue = Arc::new(ConcurrentQueue::bounded( insertion_block_size));
            let sketching_start_time = SystemTime::now();
            let sketching_start_cpu = ThreadTime::now();
            // we can create a new thread for at least nb_bases_thread_threshold bases.
            let nb_bases_thread_threshold : usize = 10_000_000;
            log::info!("threshold number of bases for thread cretaion : {:?}", nb_bases_thread_threshold);
            // a bounded blocking queue to limit the number of threads to pool_nb_thread.
            // at thread creation we send a msg into queue, at thread end we receive a msg.
            // So the length of the channel is number of active thread
            let (thread_token_sender, thread_token_receiver) = crossbeam_channel::bounded::<u32>(nb_max_threads);
            // how many msg we received
            let mut nb_msg_received = 0;
            // we must read messages, sketch and insert into hnsw
            let mut read_more = true;
            let mut nb_base_in_queue = 0;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                if !insertion_queue.is_full() {
                    let res_receive = receive.recv();
                    match res_receive {
                        Err(RecvError) =>   {   read_more = false;
                                                log::debug!("end of collector reception");
                        }
                        Ok(idsequences) => {
                            // concat the new idsketch in insertion queue.
                            nb_msg_received += 1;
                            nb_base_in_queue = idsequences.iter().fold(0, |acc, s| acc + s.get_seq_len());
                            log::debug!("nb_base_in_queue : {}", nb_base_in_queue);
                            log::debug!("read received nb seq : {}, total nb msg received : {}", idsequences.len(), nb_msg_received);
                            let p_res = insertion_queue.push(idsequences);
                            if let Err(e) = p_res {
                                match e {
                                    PushError::Full(_)   => { log::error!("queue is full"); },
                                    PushError::Closed(_) => { log::error!("queue is closed"); },
                                };
                                // 
                                log::error!("aborting in thread dispatcher, cannot push in insertion queue");
                                std::process::abort();
                            };
                        }
                    }
                }
                if read_more == false && insertion_queue.len() == 0 {
                    break;
                }
                // if insertion_queue is beyond threshold size in number of bases we can go to threaded sketching and threading insertion
                if !thread_token_sender.is_full() && (nb_base_in_queue > nb_bases_thread_threshold || insertion_queue.is_full()) {
                    let collect_sender_clone = collect_sender.clone();
                    let sketcher_clone = sketcher.clone();
                    let mut local_queue = Vec::<Vec<IdSeq> >::with_capacity(insertion_queue.len());
                    let q_len = insertion_queue.len();
                    log::debug!("insertion_queue.len() : {}", q_len);
                    for _ in 0..q_len {
                        let res_pop = insertion_queue.pop();
                        if res_pop.is_ok() {
                            local_queue.push(res_pop.unwrap());
                        }
                    }
                    nb_base_in_queue = 0; 
                    assert!(local_queue.len() > 0);
                    let _ = thread_token_sender.send(1);
                    log::debug!("nb sketching threads running = {:?}", thread_token_sender.len());
                    let thread_token_receiver_cloned = thread_token_receiver.clone();
                    scope.spawn(  move |_|  {
                        log::trace!("spawning thread on nb files : {}", local_queue.len());
                        match block_processing {
                            true => {
                                for v in &local_queue {
                                    assert_eq!(v.len(), 1);
                                }
                                let sequencegroup_ref : Vec<&SequenceAA> = local_queue.iter().map(|v| v[0].get_sequence_aa().unwrap()).collect();
                                log::debug!("calling sketch_compressedkmer nb seq : {}", sequencegroup_ref.len());
                                // computes hash signature , as we treat the file globally, the signature is indexed by filerank!
                                let signatures = sketcher_clone.sketch_compressedkmeraa(&sequencegroup_ref, kmer_hash_fn);
                                let seq_rank :  Vec<usize> = local_queue.iter().map(|v| v[0].get_filerank()).collect();
                                assert_eq!(signatures.len(), seq_rank.len(), "signatures len != seq rank len");
                                for i in 0..signatures.len() {
                                    let item: ItemDict = ItemDict::new(Id::new(local_queue[i][0].get_path(), local_queue[i][0].get_fasta_id()), local_queue[i][0].get_seq_len());
                                    let msg = CollectMsg::new((signatures[i].clone(), seq_rank[i]), item);
                                    let _ = collect_sender_clone.send(msg);
                                }
                            }
                            false => { // means we are in seq by seq mode inside a file. We get one signature by msg (or file)
                                // TODO: This can be further // either by iter or explicit thread
                                for i in 0..local_queue.len() {
                                    let sequencegroup_ref : Vec<&SequenceAA> = local_queue[i].iter().map(|v| v.get_sequence_aa().unwrap()).collect();
                                    let seq_len = sequencegroup_ref.iter().fold(0, | acc, s | acc + s.size());
                                    // as we treat the file globally, the signature is indexed by filerank!
                                    let seq_rank :  Vec<usize> = local_queue[i].iter().map(|v| v.get_filerank()).collect();
                                    let signatures = sketcher_clone.sketch_compressedkmeraa_seqs(&sequencegroup_ref, kmer_hash_fn);
                                    log::debug!("msg num : {}, nb sub seq : {}, path : {:?}, file rank : {}", i, local_queue[i].len(), local_queue[i][0].get_path(), seq_rank[0]);
                                    assert_eq!(signatures.len(), 1);
                                    // we get the item for the first seq (all sub sequences have same identity)
                                    let item: ItemDict = ItemDict::new(Id::new(local_queue[i][0].get_path(), local_queue[i][0].get_fasta_id()), seq_len);
                                    let msg = CollectMsg::new((signatures[0].clone(), seq_rank[0]), item);
                                    let _ = collect_sender_clone.send(msg);
                                }
                            }
                        }
                        // cleaning
                        local_queue.clear();
                        drop(local_queue); // a precaution as  local_queue has been moved into the thread 
                        // we free a token in queue
                        let res = thread_token_receiver_cloned.recv();
                        log::debug!("thread_token_receiver_cloned.len = {:?}", thread_token_receiver_cloned.len());
                        match res {
                            Ok(_)   => { log::debug!("thread ending OK");},
                            Err(_)  => { log::error!("thread exit with a token error, aborting");
                                        std::process::abort();
                                    },
                        }
                    });   // end internal thread 
                }
            } // end while
            // information on how time is spent
            let thread_cpu_time = sketching_start_cpu.try_elapsed();
            if thread_cpu_time.is_ok() {
                log::info!("sketcher processed in  system time(s) : {}, cpu time(s) : {}", 
                    sketching_start_time.elapsed().unwrap().as_secs(), thread_cpu_time.unwrap().as_secs());
            }
            //
            drop(receive);
            drop(collect_sender);
        }); // end of receptor thread

        // a collector task to synchronize access to hnsw and SeqDict
        scope.spawn(|_| {
            let mut msg_store = Vec::<(VecSig<Sketcher,Kmer>, usize)>::with_capacity(3 * insertion_block_size);
            let mut itemv =  Vec::<ItemDict>::with_capacity(3 * insertion_block_size);
            let mut dict_size = seqdict.get_nb_entries();
            let mut read_more = true;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = collect_receiver.recv();
                match res_receive {
                    Err(RecvError) =>   {   read_more = false;
                                            log::debug!("end of collector reception");
                    }
                    Ok(to_insert) => {
                        msg_store.push(to_insert.skecth_and_rank);
                        itemv.push(to_insert.item);
                        nb_received += 1;
                        log::debug!("collector received nb_received : {}", nb_received);
                    }
                }
                if read_more == false && msg_store.len() == 0 {
                    break;
                }
                if read_more == false || msg_store.len() > insertion_block_size {
                    log::debug!("inserting block in hnsw, nb new points : {:?}", msg_store.len());
                    let mut data_for_hnsw = Vec::<(&VecSig<Sketcher,Kmer>, usize)>::with_capacity(msg_store.len());
                    for i in 0..msg_store.len() {
                        log::debug!("inserting data id(filerank) : {}, itemdict : {}", msg_store[i].1, itemv[i].get_id().get_path());
                        // due to threading file arrive in random order, so it this thread responsability to affect id in dictionary
                        // otherwise we must introduce an indexmap in SeqDict
                        data_for_hnsw.push((&msg_store[i].0, dict_size));
                        dict_size += 1;
                    }
                    hnsw.parallel_insert(&data_for_hnsw);  
                    seqdict.0.append(&mut itemv); 
                    assert_eq!( seqdict.get_nb_entries(), hnsw.get_nb_point());
                    log::debug!(" dictionary size : {} hnsw nb points : {}", seqdict.get_nb_entries(), hnsw.get_nb_point());
                    msg_store.clear();
                    itemv.clear();
                }
            }
            //
            log::debug!("collector thread dumping hnsw , received nb_received : {}", nb_received);
            let _ = dumpall(dump_path_ref, &hnsw, &seqdict, &processing_params);
        }); // end of collector thread

    });  // end of pool
    //
    log::info!("sketchandstore, nb_sent = {}, nb_received = {}", nb_sent, nb_received);
    if nb_sent != nb_received {
        log::warn!("an error occurred  nb msg sent : {}, nb msg received : {}", nb_sent, nb_received);
    }
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


pub fn aa_process_tohnsw(dirpath : &PathBuf, filter_params : &FilterParams, processing_parameters : &ProcessingParams, others_params : &ComputingParams) {
    //
    let sketch_params = processing_parameters.get_sketching_params();
    let kmer_size = sketch_params.get_kmer_size();
    //
    match sketch_params.get_algo() {
        SketchAlgo::PROB3A => {
            //
            if kmer_size <= 6 {
                let sketcher = ProbHash3aSketch::<KmerAA32bit>::new(sketch_params);
                sketchandstore_dir_compressedkmer::<KmerAA32bit, ProbHash3aSketch::<KmerAA32bit>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size <= 12 {
                let sketcher = ProbHash3aSketch::<KmerAA64bit>::new(sketch_params);
                sketchandstore_dir_compressedkmer::<KmerAA64bit, ProbHash3aSketch::<KmerAA64bit>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            } else  {
                panic!("kmer for Amino Acids must be less or equal to 12");
            }
        }
        SketchAlgo::SUPER => {
            if kmer_size <= 6 {
                let sketcher = SuperHashSketch::<KmerAA32bit, f32>::new(sketch_params);
                sketchandstore_dir_compressedkmer::<KmerAA32bit, SuperHashSketch::<KmerAA32bit, f32>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size <= 12 {
                let sketcher = SuperHashSketch::<KmerAA64bit, f32>::new(sketch_params);
                sketchandstore_dir_compressedkmer::<KmerAA64bit, SuperHashSketch::<KmerAA64bit, f32>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);                
            }
            else {
                panic!("kmer for Amino Acids must be less or equal to 12");
            }
        }
        SketchAlgo::HLL => {
            let mut hll_params = SetSketchParams::default();
            if hll_params.get_m() < sketch_params.get_sketch_size() as u64 {
                log::warn!("!!!!!!!!!!!!!!!!!!!  need to adjust hll parameters!");
            }
            hll_params.set_m(sketch_params.get_sketch_size());
            //
            if kmer_size <= 6 {
                let sketcher = HyperLogLogSketch::<KmerAA32bit, u16>::new(sketch_params, hll_params);
                sketchandstore_dir_compressedkmer::<KmerAA32bit, HyperLogLogSketch::<KmerAA32bit, u16>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size <= 12 {
                let sketcher =  HyperLogLogSketch::<KmerAA64bit, u16>::new(sketch_params, hll_params);
                sketchandstore_dir_compressedkmer::<KmerAA64bit, HyperLogLogSketch::<KmerAA64bit, u16>>(&dirpath, sketcher, &filter_params, &processing_parameters, others_params);                
            }
            else {
                panic!("kmer for Amino Acids must be less or equal to 12");
            }
        }
        SketchAlgo::SUPER2 => {
            panic!("SUPER2 not yet implemented over AA sketching");
        }
    }
} // end of aa_process_tohnsw
