//#![allow(dead_code)]
//#![allow(unused_variables)]

//! tohnsw --dir [-d] dir --sketch [-s] size --nbng [-n] nb
//! 
//! --dir : the name of directory containing tree of GCF and GCA files 
//! --sketch gives the size of probminhash sketch ()integer value)
//! --kmer [-k] gives the size of kmer to use for generating probminhash (integer value)
//! --nbng [-n] gives the number of neihbours required in hnsw construction
// must loop on sub directories , open gzipped files
// extracts complete genomes possiby many in one file (get rid of capsid records if any)
// compute probminhash sketch and store in a Hnsw.

// one thread should read sequences and do the probminhash
// another process should store in hnsw

// hnsw should also run in a query server mode after insertion.

use clap::{App, Arg};

use std::path::Path;

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

//
use std::time::{SystemTime};
use cpu_time::ProcessTime;

// for multithreading
use crossbeam_channel::*;

// our crate
use hnsw_rs::prelude::*;
use kmerutils::base::{sequence::*, Kmer32bit, KmerT};
use kmerutils::sketching::*;
use kmerutils::sketching::seqsketchjaccard::SeqSketcher;


use archaea::utils::idsketch::{SeqDict,Id, IdSeq};

use archaea::utils::files::{process_dir,process_file};




struct HnswParams {
    #[allow(dead_code)]
    nbng : usize,
    /// expected number of sequences to store
    capacity : usize,
    ef_search : usize,
    max_nb_conn : usize,
}

// install a logger facility
pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}

// this function does the sketching and hnsw store of a whole directory
fn sketchandstore_dir(dirpath : &Path, sketcher_params : &SeqSketcher, hnsw_params : &HnswParams) {
    //
    log::trace!("sketchandstore_dir processing dir {}", dirpath.to_str().unwrap());
    log::info!("sketchandstore_dir {}", dirpath.to_str().unwrap());
    let start_t = SystemTime::now();
    let cpu_start = ProcessTime::now();
    //
    // a queue of signature waiting to be inserted , size must be sufficient to benefit from threaded probminhash and insert
    let insertion_block_size = 5000;
    let mut insertion_queue : Vec<IdSeq>= Vec::with_capacity(insertion_block_size);
    // TODO must get ef_search from clap via hnswparams
    let hnsw = Hnsw::<u32, DistHamming>::new(hnsw_params.max_nb_conn , hnsw_params.capacity , 16, hnsw_params.ef_search, DistHamming{});
    //
    // Sketcher allocation, we do not need reverse complement hashing as we sketch assembled genomes. (Jianshu Zhao)
    //
    let kmer_hash_fn = | kmer : &Kmer32bit | -> u32 {
        let value_mask :u32 = (0b1 << (2*kmer.get_nb_base())) - 1;
        let hashval = kmer.0 & value_mask;
        hashval
    };
    let sketcher = seqsketchjaccard::SeqSketcher::new(sketcher_params.get_kmer_size(), sketcher_params.get_sketch_size());
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(10_000);
    // launch process_dir in a thread or async
    crossbeam_utils::thread::scope(|scope| {
        // sequence sending, productor thread
        let mut nb_sent = 0;
        let sender_handle = scope.spawn(move |_|   {
            let res_nb_sent = process_dir(dirpath, &process_file, &send);
            match res_nb_sent {
                Ok(nb_really_sent) => {
                    nb_sent = nb_really_sent;
                    println!("process_dir processed nb sequences : {}", nb_sent);
                }
                Err(_) => {
                    println!("some error occured in process_dir");
                }
            };
            drop(send);
            Box::new(nb_sent)
        });
        // sequence reception, consumer thread
        let receptor_handle = scope.spawn(move |_| {
            let mut seqdict = SeqDict::new(100000);
            // we must read messages, sketch and insert into hnsw
            let mut read_more = true;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = receive.recv();
                match res_receive {
                    Err(RecvError) => { read_more = false;
                        // sketch the content of  insertion_queue if not empty
                        if insertion_queue.len() > 0 {
                            let sequencegroup_ref : Vec<&Sequence> = insertion_queue.iter().map(|s| s.get_sequence()).collect();
                            // collect rank
                            let seq_rank :  Vec<usize> = insertion_queue.iter().map(|s| s.get_rank()).collect();
                            // collect Id
                            let mut seq_id :  Vec<Id> = insertion_queue.iter().map(|s| Id::new(s.get_path(), s.get_fasta_id())).collect();
                            seqdict.0.append(&mut seq_id);
                            let signatures = sketcher.sketch_probminhash3a_kmer32bit(&sequencegroup_ref, kmer_hash_fn);                            
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch for hnsw insertion
                            let mut data_for_hnsw = Vec::<(&Vec<u32>, usize)>::with_capacity(signatures.len());
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
                            let sequencegroup_ref : Vec<&Sequence> = insertion_queue.iter().map(|s| s.get_sequence()).collect();
                            let seq_rank :  Vec<usize> = insertion_queue.iter().map(|s| s.get_rank()).collect();
                            // collect Id
                            let mut seq_id :  Vec<Id> = insertion_queue.iter().map(|s| Id::new(s.get_path(), s.get_fasta_id())).collect();
                            seqdict.0.append(&mut seq_id);
                            // computes hash signature
                            let signatures = sketcher.sketch_probminhash3a_kmer32bit(&sequencegroup_ref, kmer_hash_fn);
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch, inserting unique id, for hnsw insertion
                            let mut data_for_hnsw = Vec::<(&Vec<u32>, usize)>::with_capacity(signatures.len());
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
            // and finally dump sketchparams
            let _ = sketcher_params.dump_json(&"sketchparams_dump.json".to_string());
            //
            Box::new(seqdict.0.len())
        }); // end of receptor thread
        // now we must join handles
        let nb_sent = sender_handle.join().unwrap();
        let nb_received = receptor_handle.join().unwrap();
        log::debug!("sketchandstore, nb_sent = {}, nb_received = {}", nb_sent, nb_received);
        if nb_sent != nb_received {
            log::error!("an error occurred  nb msg sent : {}, nb msg received : {}", nb_sent, nb_received);
        }
    }).unwrap();  // end of scope
    //
    let cpu_time = cpu_start.elapsed().as_secs();
    log::info!("process_dir : cpu time(s) {}", cpu_time);
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
    log::info!("process_dir : elapsed time(s) {}", elapsed_t);
} // end of sketchandstore




fn main() {
    let _ = init_log();
    //
    let matches = App::new("tohnsw")
        .arg(Arg::with_name("dir")
            .long("dir")
            .short("d")
            .takes_value(true)
            .help("name of directory containing genomes to index"))
        .arg(Arg::with_name("kmer_size")
            .long("kmer")
            .short("k")
            .takes_value(true)
            .help("expecting a kmer size"))
        .arg(Arg::with_name("sketch size")
            .long("sketch")
            .short("s")
            .default_value("8")
            .help("size of probinhash sketch, default to 8"))
        .arg(Arg::with_name("neighbours")
            .long("nbng")
            .short("n")
            .takes_value(true)
            .help("must specify number of neighbours in hnsw"))
        .get_matches();

    // decode matches, check for dir
        let datadir;
        if matches.is_present("dir") {
            println!("decoding argument dir");
            datadir = matches.value_of("dir").ok_or("").unwrap().parse::<String>().unwrap();
            if datadir == "" {
                println!("parsing of dir failed");
                std::process::exit(1);
            }
        }
        else {
            println!("-d dirname is mandatory");
            std::process::exit(1);
        }
        let dirpath = Path::new(&datadir);
        if !dirpath.is_dir() {
            println!("error not a directory : {:?}", datadir);
            std::process::exit(1);
        }
        // get sketching params
        let mut sketch_size = 256;
        if matches.is_present("size") {
            sketch_size = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("sketching size {}", sketch_size);
        }
        else {
            println!("using default sketch size {}", sketch_size);
        }
        //
        let mut kmer_size = 6;
        if matches.is_present("kmer_size") {
            kmer_size = matches.value_of("kmer_size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("kmer size {}", kmer_size);
        }
        else {
            println!("using default kmer size {}", kmer_size);
        }
        let sketch_params =  SeqSketcher::new(kmer_size as usize, sketch_size as usize);  
        //
        let nbng;
        if matches.is_present("neighbours") {
            nbng = matches.value_of("neighbours").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("nb neighbours you will need in hnsw requests {}", nbng);
        }        
        else {
            println!("-n nbng is mandatory");
            std::process::exit(1);
        }
        // max_nb_conn must be adapted to the number of neighbours we will want in searches.
        let max_nb_conn = 128.min(2 * nbng as usize);
        let ef_search = 200;
        let hnswparams = HnswParams{nbng : nbng as usize, capacity : 700000, ef_search, max_nb_conn};
        //
        //
        sketchandstore_dir(&dirpath, &sketch_params, &hnswparams);


        //
 } // end of main