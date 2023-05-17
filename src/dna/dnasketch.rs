//! dna sketching
//! 
//! 
//! 

use std::time::{SystemTime};
use cpu_time::{ProcessTime,ThreadTime};


use std::path::{PathBuf};

// for multithreading
use std::sync::Arc;
use crossbeam_channel::*;
use concurrent_queue::{ConcurrentQueue, PushError};


use serde::{de::DeserializeOwned, Serialize};

use std::fmt::{Debug};

// our crate
use hnsw_rs::prelude::*;


use kmerutils::base::{kmergenerator::*, Kmer32bit, Kmer64bit, CompressedKmerT};
use kmerutils::sketching::seqsketchjaccard::{SeqSketcherT, ProbHash3aSketch, SuperHashSketch, HyperLogLogSketch, SuperHash2Sketch};
use kmerutils::sketcharg::DataType;

use probminhash::{setsketcher::SetSketchParams};

use crate::utils::{idsketch::*, reloadhnsw};
use crate::utils::files::{process_dir,process_dir_parallel, ProcessingState};

use crate::dna::dnafiles::{process_file_in_one_block, process_file_by_sequence, process_buffer_in_one_block, process_buffer_by_sequence};

use crate::utils::parameters::*;


// Sig is basic item of a signature , VecSig is a vector of such items
type VecSig<Sketcher, Kmer>  = Vec< <Sketcher as SeqSketcherT<Kmer>>::Sig>;


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

// hnsw_pb  must contains directory of hnsw database.
// In creation mode it is the directory containings file to proceess. In add mode it is where hnsw database reside,
// and the directory containing new data are taken from computingParams if in add mode
fn sketchandstore_dir_compressedkmer<Kmer:CompressedKmerT+KmerBuilder<Kmer>, Sketcher : SeqSketcherT<Kmer>>(hnsw_pb : &PathBuf, sketcher : Sketcher ,
                    filter_params: &FilterParams, 
                    processing_params : &ProcessingParams, 
                    other_params : &ComputingParams) 

        where   Sketcher : SeqSketcherT<Kmer> + Clone + Send + Sync,
                <Sketcher as SeqSketcherT<Kmer>>::Sig : 'static + Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
                KmerGenerator<Kmer> :  KmerGenerationPattern<Kmer>, 
                DistHamming : Distance< <Sketcher as  SeqSketcherT<Kmer>>::Sig>,
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
        true => { 5000.min(2 * other_params.get_nb_files_par()) },
        _    => { 5000 },
    };
    log::debug!("insertion_block_size : {}", insertion_block_size);
    //
    let mut hnsw : Hnsw::< <Sketcher as SeqSketcherT<Kmer>>::Sig, DistHamming>;
    let mut state :ProcessingState;
    //
    let toprocess_path : PathBuf;
    //
    if other_params.get_adding_mode() {
        let toprocess_str = other_params.get_add_dir();
        toprocess_path = PathBuf::from(toprocess_str);
        // in this case we must reload
        log::info!("dnasketch::sketchandstore_dir_compressedkmer will add new data from from directory {:?} to hnsw dir {:?}", hnsw_pb, toprocess_path);
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
        hnsw = Hnsw::< <Sketcher as SeqSketcherT<Kmer>>::Sig, DistHamming>::new(hnsw_params.get_max_nb_connection() as usize , hnsw_params.capacity , 16, hnsw_params.get_ef(), DistHamming{});
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
            crossbeam_channel::bounded::<CollectMsg<<Sketcher as SeqSketcherT<Kmer>>::Sig>>(insertion_block_size+1);
    //
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let pool_nb_thread = pool.current_num_threads();
    log::info!("nb threads in pool : {:?}", pool_nb_thread);

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
                                        nb_files_by_group, &process_buffer_in_one_block, &send);
                    } // end case parallel io
                    else {
                        log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir serial");
                        res_nb_sent = process_dir(&mut state, &DataType::DNA, &toprocess_path, filter_params, 
                                        &process_file_in_one_block, &send);
                    }                
                }
                false => {
                    log::info!("dnasketch::sketchandstore_dir_compressedkmer : seq by seq processing");
                    if other_params.get_parallel_io() {
                        let nb_files_by_group = other_params.get_nb_files_par();
                        log::info!("dnasketch::sketchandstore_dir_compressedkmer : calling process_dir_parallel, nb_files in parallel : {}", nb_files_by_group);
                        res_nb_sent = process_dir_parallel(&mut state, &DataType::DNA,  &toprocess_path, filter_params, 
                                nb_files_by_group, &process_buffer_by_sequence, &send);
                    }
                    else {
                        log::info!("processing by sequence, sketching whole file globally");
                        res_nb_sent = process_dir(&mut state, &DataType::DNA, &toprocess_path, filter_params,
                                     &process_file_by_sequence, &send);
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
            let insertion_queue = Arc::new(ConcurrentQueue::bounded(2 * insertion_block_size));
            let sketching_start_time = SystemTime::now();
            let sketching_start_cpu = ThreadTime::now();
            // we can create a new thread for at least nb_bases_thread_threshold bases.
            let nb_bases_thread_threshold : usize = 1_000_000;
            log::info!("threshold number of bases for thread cretaion : {:?}", nb_bases_thread_threshold);
            // a bounded blocking queue to limit the number of threads to pool_nb_thread.
            // at thread creation we send a msg into queue, at thread end we receive a msg.
            // So the length of the channel is number of active thread
            let (thread_token_sender, thread_token_receiver) = crossbeam_channel::bounded::<u32>(pool_nb_thread);
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
                if nb_base_in_queue > nb_bases_thread_threshold || insertion_queue.is_full() {
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
                                let sequencegroup_ref : Vec<&Sequence> = local_queue.iter().map(|v| v[0].get_sequence_dna().unwrap()).collect();
                                log::debug!("calling sketch_compressedkmer nb seq : {}", sequencegroup_ref.len());
                                // computes hash signature , as we treat the file globally, the signature is indexed by filerank!
                                let signatures = sketcher_clone.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
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
                                    let sequencegroup_ref : Vec<&Sequence> = local_queue[i].iter().map(|v| v.get_sequence_dna().unwrap()).collect();
                                    let seq_len = sequencegroup_ref.iter().fold(0, | acc, s | acc + s.size());
                                    // as we treat the file globally, the signature is indexed by filerank!
                                    let seq_rank :  Vec<usize> = local_queue[i].iter().map(|v| v.get_filerank()).collect();
                                    let signatures = sketcher_clone.sketch_compressedkmer_seqs(&sequencegroup_ref, kmer_hash_fn);
                                    log::debug!("msg num : {}, nb sub seq : {}, path : {:?}", i, local_queue[i].len(), local_queue[i][0].get_path());
                                    assert_eq!(signatures.len(), 1);
                                    // we get the item for the first seq (all sub sequences have same identity)
                                    let item: ItemDict = ItemDict::new(Id::new(local_queue[i][0].get_path(), local_queue[i][0].get_fasta_id()), seq_len);
                                    let msg = CollectMsg::new((signatures[i].clone(), seq_rank[i]), item);
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
                if read_more == false || msg_store.len() > insertion_block_size {
                    log::debug!("inserting block in hnsw, nb new points : {:?}", msg_store.len());
                    let mut data_for_hnsw = Vec::<(&VecSig<Sketcher,Kmer>, usize)>::with_capacity(msg_store.len());
                    for i in 0..msg_store.len() {
                        data_for_hnsw.push((&msg_store[i].0, msg_store[i].1));
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
/// dirpath is where hnsw base will be
pub fn dna_process_tohnsw(dirpath : &PathBuf, filter_params : &FilterParams, processing_parameters : &ProcessingParams, others_params : &ComputingParams) {
    // dispatch according to kmer_size
    let kmer_size = processing_parameters.get_sketching_params().get_kmer_size();
    //
    let sketchparams = processing_parameters.get_sketching_params();
    match sketchparams.get_algo() {
        SketchAlgo::PROB3A => {
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = ProbHash3aSketch::<Kmer32bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer32bit, ProbHash3aSketch::<Kmer32bit> >(dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = ProbHash3aSketch::<Kmer16b32bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, ProbHash3aSketch::<Kmer16b32bit>>(dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = ProbHash3aSketch::<Kmer64bit>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer64bit, ProbHash3aSketch::<Kmer64bit>>(dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        } // end PPROB
        SketchAlgo::SUPER => {
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = SuperHashSketch::<Kmer32bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer32bit, SuperHashSketch::<Kmer32bit, f32> >(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = SuperHashSketch::<Kmer16b32bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, SuperHashSketch::<Kmer16b32bit, f32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = SuperHashSketch::<Kmer64bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer64bit, SuperHashSketch::<Kmer64bit, f32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        }
        SketchAlgo::HLL => {
            // TODO !! adjust hll parameters to 
            let mut hll_params = SetSketchParams::default();
            if hll_params.get_m() < sketchparams.get_sketch_size() as u64 {
                log::warn!("!!!!!!!!!!!!!!!!!!!  need to adjust hll parameters!");
            }
            hll_params.set_m(sketchparams.get_sketch_size());
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = HyperLogLogSketch::<Kmer32bit, u16>::new(sketchparams, hll_params);
                sketchandstore_dir_compressedkmer::<Kmer32bit, HyperLogLogSketch::<Kmer32bit, u16> >(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = HyperLogLogSketch::<Kmer16b32bit, u16>::new(sketchparams, hll_params);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, HyperLogLogSketch::<Kmer16b32bit, u16>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = HyperLogLogSketch::<Kmer64bit, u16>::new(sketchparams, hll_params);
                sketchandstore_dir_compressedkmer::<Kmer64bit, HyperLogLogSketch::<Kmer64bit, u16>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        } // end match HLL
        SketchAlgo::SUPER2 => {
            if kmer_size <= 14 {
                // allocated the correct sketcher. SuperHash2Sketch should enable BuildHasher choice to sketch to u32.
                // we used just in tests
                let bh = std::hash::BuildHasherDefault::<fxhash::FxHasher32>::default();
                let sketcher = SuperHash2Sketch::<Kmer32bit, u32, fxhash::FxHasher32>::new(sketchparams, bh);
                sketchandstore_dir_compressedkmer::<Kmer32bit, SuperHash2Sketch::<Kmer32bit, u32, fxhash::FxHasher32> >(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let bh = std::hash::BuildHasherDefault::<fxhash::FxHasher32>::default();
                let sketcher = SuperHash2Sketch::<Kmer16b32bit, u32, fxhash::FxHasher32>::new(sketchparams, bh);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, SuperHash2Sketch::<Kmer16b32bit, u32, fxhash::FxHasher32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let bh = std::hash::BuildHasherDefault::<fxhash::FxHasher64>::default();
                let sketcher = SuperHash2Sketch::<Kmer64bit, u64, fxhash::FxHasher64>::new(sketchparams, bh);
                sketchandstore_dir_compressedkmer::<Kmer64bit, SuperHash2Sketch::<Kmer64bit, u64, fxhash::FxHasher64>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        }        
    }
} // end of dna_process

//===========================================================================
// TODO to be transferred into utils and factored for dna and aa
// Dump reload stuff

// This function dumps hnsw , seqdict and processing params in the same directory given by dump_path_ref
fn dumpall<Sig>(dump_path_ref : &PathBuf, hnsw : &Hnsw<Sig, DistHamming>, seqdict : &SeqDict, processing_params : &ProcessingParams) -> anyhow::Result<(),String>
    where Sig : Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug ,
    DistHamming : Distance<Sig> {

    if  hnsw.get_nb_point() > 0 {
        let mut hnsw_dump = dump_path_ref.to_path_buf().clone();
        hnsw_dump.push("hnswdump");
        let hnswdumpname = String::from(hnsw_dump.to_str().unwrap());
        log::info!("going to dump hnsw with prefix : {:?}", hnswdumpname);
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
        let mut seq_pb = dump_path_ref.clone();
        seq_pb.push("seqdict.json");
        let seqdict_name = String::from(seq_pb.to_str().unwrap());
        let resdump = seqdict.dump(seqdict_name);
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
    let _ = processing_params.dump_json(dump_path_ref);
    //
    Ok(())
} // end of dumpall



fn reload_seqdict(dump_path_ref : &PathBuf) -> SeqDict {
       // must reload seqdict
       let mut filepath = PathBuf::from(dump_path_ref.clone());
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
           let seqdict = res_reload.unwrap();
           return seqdict;
       }    
} // end of reload_seqdict



// retrieve or allocate a SeqDict depending on use case
fn get_seqdict(dump_path_ref : &PathBuf, other_params : &ComputingParams) -> anyhow::Result<SeqDict> {
    //
    let seqdict = if other_params.get_adding_mode() {
        reload_seqdict(dump_path_ref)
    }
    else {
        SeqDict::new(100000)
    };
    //
    return Ok(seqdict);
} // end of get_seqdict