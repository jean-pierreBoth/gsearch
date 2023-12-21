//! dna sketching
//! 
//! 
//! 

use std::time::SystemTime;
use cpu_time::{ProcessTime,ThreadTime};
use std::time::Duration;

use std::path::PathBuf;

// for multithreading
use std::sync::Arc;
use crossbeam_channel::*;
use crossbeam::sync::Parker;
use concurrent_queue::{ConcurrentQueue, PushError};


use serde::{de::DeserializeOwned, Serialize};

use std::fmt::Debug;

// our crate
use hnsw_rs::{prelude::*, hnswio::HnswIo};


use kmerutils::base::{kmergenerator::*, Kmer32bit, Kmer64bit, CompressedKmerT};
use kmerutils::sketching::setsketchert::*;
use kmerutils::sketcharg::DataType;

use probminhash::setsketcher::SetSketchParams;

use crate::utils::{idsketch::*, reloadhnsw};
use crate::utils::dumpload::*;
use crate::utils::files::{process_dir,process_dir_parallel, ProcessingState};

use crate::dna::dnafiles::{process_file_in_one_block, process_file_by_sequence, process_buffer_in_one_block, process_buffer_by_sequence};

use crate::utils::parameters::*;



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
                    computing_params : &ComputingParams) 

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
    // processing parameters, blocked or by sequence , threads , parallel io...
    // a queue of signature waiting to be inserted , size must be sufficient to benefit from threaded probminhash and insert
    // and not too large to spare memory. If parallel_io is set dimension message queue to size of group
    // for files of size more than Gb we must use pario to limit memory, but leave enough msg in queue to get // sketch and insertion 
    // The number of threads must be adapted to memory size when sketching files of size tens of giabytes
    let block_processing = processing_params.get_block_flag();
    let insertion_block_size = match computing_params.get_parallel_io() {
        true => { 5000.min(computing_params.get_nb_files_par()) },
        _    => { 2000 },
    };
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let pool_nb_thread = pool.current_num_threads();
    let mut nb_max_threads = computing_params.get_sketching_nbthread();
    if nb_max_threads == 0 {
        nb_max_threads = insertion_block_size.min(pool_nb_thread).max(1);
    }
    log::debug!("insertion_block_size : {}", insertion_block_size);
    log::info!("nb threads in pool : {:?}, using nb threads : {}", pool_nb_thread, nb_max_threads);
    //
    let mut state :ProcessingState;
    //
    let toprocess_path : PathBuf;
    //
    let mut hnswio : HnswIo;
    let mut hnsw : Hnsw::< <Sketcher as SeqSketcherT<Kmer>>::Sig, DistHamming>;
    //
    if computing_params.get_adding_mode() {
        let toprocess_str = computing_params.get_add_dir();
        toprocess_path = PathBuf::from(toprocess_str);
        // in this case we must reload
        log::info!("dnasketch::sketchandstore_dir_compressedkmer will add new data from directory {:?} to hnsw dir {:?}", toprocess_path, hnsw_pb);

        let hnswio_res = reloadhnsw::get_hnswio(hnsw_pb);
        if hnswio_res.is_err() {
            std::panic!("cannot reload from directory {:?} failed msg : {:?}", &hnsw_pb, hnswio_res.err());
        }
        hnswio = hnswio_res.unwrap();
        let hnsw_res = hnswio.load_hnsw();
        if hnsw_res.is_err() {
            log::error!("cannot reload hnsw from directory : {:?}", &hnsw_pb);
            std::process::exit(1);
        }
        else { 
            hnsw = hnsw_res.unwrap();
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
    seqdict = get_seqdict(&hnsw_pb, computing_params).unwrap();
    //
    // where do we dump hnsw* seqdict and so on
    // If in add mode we dump where is already an hnsw database
    // If creation mode we dump in .
    //
    let dump_path= if computing_params.get_adding_mode() {
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
    // capacity not too large beccause of memory. If data is large cpu takes more time
    // than reading. The sketching queue is more important 
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(1.max(insertion_block_size));
    // to send sketch result to a collector task
    let (collect_sender , collect_receiver) = 
            crossbeam_channel::bounded::<CollectMsg<<Sketcher as SeqSketcherT<Kmer>>::Sig>>(insertion_block_size+1);
    //

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
                    if computing_params.get_parallel_io() {
                        let nb_files_by_group = computing_params.get_nb_files_par();
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
                    if computing_params.get_parallel_io() {
                        let nb_files_by_group = computing_params.get_nb_files_par();
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
            let sketcher_queue = Arc::new(ConcurrentQueue::bounded( insertion_block_size));
            let sketching_start_time = SystemTime::now();
            let sketching_start_cpu = ThreadTime::now();
            let parker: Parker = Parker::new();
            // we can create a new thread for at least nb_bases_thread_threshold bases.
            let nb_bases_thread_threshold : usize = 100_000_000;
            log::info!("threshold number of bases for thread creation : {:?}", nb_bases_thread_threshold);
            // a bounded blocking queue to limit the number of threads to nb_max_threads
            // at thread creation we send a msg into queue, at thread end we receive a msg.
            // So the length of the channel is number of active thread
            let (thread_token_sender, thread_token_receiver) = crossbeam_channel::bounded::<u32>(nb_max_threads);
            // how many msg we received
            let mut nb_msg_received = 0;
            // we must read messages, sketch and insert into hnsw
            let more = true;
            let mut read_more = true;
            let mut nb_base_in_queue= 0;
            while more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                if !sketcher_queue.is_full() {
                    let res_receive = receive.recv();
                    match res_receive {
                        Err(RecvError) =>   {  
                                if read_more {
                                    // first time we reach that state we log it
                                    log::info!("end of file reception");
                                }
                                read_more = false;
                        }
                        Ok(idsequences) => {
                            // concat the new idsketch in insertion queue.
                            nb_msg_received += 1;
                            nb_base_in_queue += idsequences.iter().fold(0, |acc, s| acc + s.get_seq_len());
                            log::debug!("nb_base_in_queue : {}", nb_base_in_queue);
                            log::debug!("read received nb seq : {}, total nb msg received : {}", idsequences.len(), nb_msg_received);
                            let p_res = sketcher_queue.push(idsequences);
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
                // 
                if read_more == false && sketcher_queue.len() == 0 && thread_token_sender.len() == 0 {
                    log::info!("no more read to come, no more data to sketch, no more sketcher running ...");
                    break;
                }
                else if read_more == false && sketcher_queue.len() == 0 {
                    // TODO here we just have to wait all thread end, should use a condvar, now we check every second
                    parker.park_timeout(Duration::from_millis(1000));
                }
                // if sketching_queue is beyond threshold size in number of bases we can go to threaded sketching and threading insertion
                if nb_base_in_queue >= nb_bases_thread_threshold || sketcher_queue.is_full() || (!read_more && !sketcher_queue.is_empty()) {
                    let nb_free_thread = thread_token_sender.capacity().unwrap() as i32 - thread_token_sender.len() as i32;
                    log::debug!("nb_free_thread : {}",nb_free_thread);
                    log::debug!("thread_token_sender.len() = {}, thread_token_receiver.len() : {}", thread_token_sender.len(), thread_token_receiver.len());
                    let collect_sender_clone = collect_sender.clone();
                    let sketcher_clone = sketcher.clone();
                    let mut local_queue = Vec::<Vec<IdSeq> >::with_capacity(sketcher_queue.len());
                    let to_pop = sketcher_queue.len();
                    log::debug!("popping to local_queue {}, sketching_queue.len() : {}", to_pop, sketcher_queue.len());
                    for _ in 0..to_pop {
                        let mut nb_popped_bases : usize = 0;
                        let res_pop = sketcher_queue.pop();
                        if res_pop.is_ok() {
                            nb_popped_bases += res_pop.as_ref().unwrap().iter().fold(0, |acc, s| acc + s.get_seq_len());
                            local_queue.push(res_pop.unwrap());
                        }
                        assert!(nb_base_in_queue >= nb_popped_bases);
                        nb_base_in_queue -= nb_popped_bases; 
                    }
                    log::debug!("after popping , nb_base_in_queue : {}", nb_base_in_queue);
                    assert!(local_queue.len() > 0);
                    // recall this will thread creation if too many threads already
                    let _ = thread_token_sender.send(1);
                    log::debug!("sending token sender tokens = {:?}, nb_token_reciever : {}", thread_token_sender.len(), thread_token_receiver.len());
                    let thread_token_receiver_cloned = thread_token_receiver.clone();
                    scope.spawn(  move |_|  {
                        log::trace!("spawning thread on nb files : {}", local_queue.len());
                        match block_processing {
                            true => {
                                // each msg is a concatenated seq, so unit length
                                for v in &local_queue {
                                    assert_eq!(v.len(), 1);
                                }
                                let sequencegroup_ref : Vec<&Sequence> = local_queue.iter().map(|v| v[0].get_sequence_dna().unwrap()).collect();
                                log::debug!("calling sketch_compressedkmer nb seq : {}", sequencegroup_ref.len());
                                // computes hash signature , as we treat the file globally, the signature is indexed by filerank!
                                let signatures = sketcher_clone.sketch_compressedkmer(&sequencegroup_ref, kmer_hash_fn);
                                let seq_rank :  Vec<usize> = local_queue.iter().map(|v| v[0].get_filerank()).collect();
                                assert_eq!(signatures.len(), seq_rank.len(), "signatures len != seq rank len");
                                log::debug!("sending signatures : {}", signatures.len());
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
                                    // seq_rank has size of sequencegroup_ref but all values in seq_rank should be equal to filerank
                                    let seq_rank = local_queue[i][0].get_filerank();
                                    if local_queue[i].len() > 1 {
                                        assert_eq!(seq_rank, local_queue[i][1].get_filerank());
                                    }
                                    let signatures = sketcher_clone.sketch_compressedkmer_seqs(&sequencegroup_ref, kmer_hash_fn);
                                    log::debug!("sketched msg num : {}, nb sub seq : {}, path : {:?}, file rank : {}", i, local_queue[i].len(), local_queue[i][0].get_path(), seq_rank);
                                    assert_eq!(signatures.len(), 1);
                                    // we get the item for the first seq (all sub sequences have same identity)
                                    let item: ItemDict = ItemDict::new(Id::new(local_queue[i][0].get_path(), local_queue[i][0].get_fasta_id()), seq_len);
                                    let msg = CollectMsg::new((signatures[0].clone(), seq_rank), item);
                                    let _ = collect_sender_clone.send(msg);
                                }
                            }
                        }
                        // cleaning
                        while let Some(v) = local_queue.pop() {
                            drop(v);
                        }
                        drop(local_queue); // a precaution as  local_queue has been moved into the thread 
                        drop(collect_sender_clone);
                        drop(sketcher_clone);
                        // we free a token in queue
                        let res = thread_token_receiver_cloned.recv();
                        log::debug!("thread_token_receiver_cloned.len = {:?}", thread_token_receiver_cloned.len());
                        match res {
                            Ok(_)   => { log::debug!("internal sketching thread ending OK");},
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
        }); // end of sketcher thread

        // a collector task to synchronize access to hnsw and SeqDict
        // Note that doing parallel insertion requires a new thread that put us on scheduling dependance and can cause deadlock!!!
        // so we insert at each msh received!
        scope.spawn(|_| {
            type VecSig<Sketcher, Kmer>  = Vec< <Sketcher as SeqSketcherT<Kmer>>::Sig>;
            let store_size = 100_000;
            let mut msg_store = Vec::<(VecSig<Sketcher,Kmer>, usize)>::with_capacity(store_size);
            let mut itemv =  Vec::<ItemDict>::with_capacity(store_size);
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
                        log::debug!("collector received nb_received : {}, receiver len : {}", nb_received, collect_receiver.len());
                    }
                }
                // we defer all insertion to end of sketching as we observed bad interaction between // iter and threads
                if read_more == false {
                    let insertion_start_time = SystemTime::now();
                    let insertion_cpu_start: ProcessTime = ProcessTime::now();
                    log::debug!("inserting block in hnsw, nb new points : {:?}", msg_store.len());
                    let mut data_for_hnsw = Vec::<(&VecSig<Sketcher,Kmer>, usize)>::with_capacity(msg_store.len());
                    for i in 0..msg_store.len() {
                        log::debug!("inserting data id(filerank) : {}, itemdict : {}", msg_store[i].1, itemv[i].get_id().get_path());
                        // due to threading file arrive in random order, so it this thread responsability to affect id in dictionary
                        // otherwise we must introduce an indexmap in SeqDict
                        data_for_hnsw.push((&msg_store[i].0, dict_size));
                        dict_size += 1;
                    }
                    log::info!("parallel insertion ...");
                    hnsw.parallel_insert(&data_for_hnsw);
                    log::info!("parallel insert done");  
                    seqdict.0.append(&mut itemv); 
                    assert_eq!( seqdict.get_nb_entries(), hnsw.get_nb_point());
                    // information on how time is spent
                    let thread_cpu_time = insertion_cpu_start.try_elapsed();
                    if thread_cpu_time.is_ok() {
                        log::info!("hnsw insertion processed in  system time(s) : {}, cpu time(s) : {}", 
                            insertion_start_time.elapsed().unwrap().as_secs(), thread_cpu_time.unwrap().as_secs());
                    }
                    log::debug!(" dictionary size : {} hnsw nb points : {}", seqdict.get_nb_entries(), hnsw.get_nb_point());
                    msg_store.clear();
                    itemv.clear();
                }
            }
            //
            assert_eq!(seqdict.get_nb_entries(), hnsw.get_nb_point());
            log::debug!("collector thread dumping hnsw , received nb_received : {}", nb_received);
            let _ = dumpall(dump_path_ref, &hnsw, &seqdict, &processing_params);
        }); // end of collector thread

    });  // end of pool
    //
    log::info!("sketchandstore, nb_sent = {}, nb_received = {}", nb_sent, nb_received);
    if nb_sent != nb_received {
        log::error!("==============================================================================");
        log::error!("an error occurred  nb msg sent : {}, nb msg received : {}", nb_sent, nb_received);
        log::error!("==============================================================================");
        std::process::exit(1);
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
            // we need number of threads of cpu and we have the sketching thread parameter in ComputingParams
            // if memory limits the number of sketching threads, we will use more in iterators.
            let nb_cpus = num_cpus::get();
            log::info!("nb cpus : {}", nb_cpus);
            let nb_iter_threads = ((nb_cpus.max(4) - 4) / others_params.get_sketching_nbthread().max(1)).max(1);
            let hll_seqs_threading = HllSeqsThreading::new(nb_iter_threads, 10_000_000);
            log::info!("HllSeqsThreading : {:?}", hll_seqs_threading);
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = HyperLogLogSketch::<Kmer32bit, u16>::new(sketchparams, hll_params, hll_seqs_threading);
                sketchandstore_dir_compressedkmer::<Kmer32bit, HyperLogLogSketch::<Kmer32bit, u16> >(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = HyperLogLogSketch::<Kmer16b32bit, u16>::new(sketchparams, hll_params, hll_seqs_threading);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, HyperLogLogSketch::<Kmer16b32bit, u16>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = HyperLogLogSketch::<Kmer64bit, u16>::new(sketchparams, hll_params, hll_seqs_threading);
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
        SketchAlgo::OPTDENS => {
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = OptDensHashSketch::<Kmer32bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer32bit, OptDensHashSketch::<Kmer32bit, f32> >(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = OptDensHashSketch::<Kmer16b32bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, OptDensHashSketch::<Kmer16b32bit, f32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = OptDensHashSketch::<Kmer64bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer64bit, OptDensHashSketch::<Kmer64bit, f32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        }
        //
        SketchAlgo::REVOPTDENS => {
            if kmer_size <= 14 {
                // allocated the correct sketcher
                let sketcher = RevOptDensHashSketch::<Kmer32bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer32bit, RevOptDensHashSketch::<Kmer32bit, f32> >(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if kmer_size == 16 {
                let sketcher = RevOptDensHashSketch::<Kmer16b32bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer16b32bit, RevOptDensHashSketch::<Kmer16b32bit, f32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else if  kmer_size <= 32 {
                let sketcher = RevOptDensHashSketch::<Kmer64bit, f32>::new(sketchparams);
                sketchandstore_dir_compressedkmer::<Kmer64bit, RevOptDensHashSketch::<Kmer64bit, f32>>(&dirpath, sketcher, 
                            &filter_params, &processing_parameters, others_params);
            }
            else {
                panic!("kmers cannot be greater than 32");
            }
        }      
    }
} // end of dna_process

//===========================================================================
