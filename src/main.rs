#![allow(dead_code)]
#![allow(unused_variables)]

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

use std::io;
use std::fs::{self, DirEntry};
use std::path::Path;

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

// for multithreading
use crossbeam_channel::*;

// our crate
use hnsw_rs::prelude::*;
use kmerutils::base::{sequence::*, KmerT, Kmer32bit};
use kmerutils::sketching::*;

mod idsketch;
pub use idsketch::*;


// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}


/// A compressed sequence compressed to 2 bit / base
/// This structure is used for returning info from function process_file
/// and is sent to sketcher
pub struct IdSeq {
    /// as read is sequential we can identify uniquely sequence in hnsw
    rank : usize,
    /// id of genome Sketched
    id : String,
    /// Sequence
    seq : Sequence
}  // end of IdSeq


struct SketcherParams {
    kmer_size : usize,
    sketch_size : usize,
}


struct HnswParams {
    nbng : usize,
    ef_search : usize,
    max_nb_conn : usize,
}

// returns true if file is a fasta file (possibly gzipped)
// filename are of type GCA[GCF]_000091165.1_genomic.fna.gz
fn is_fasta_file(file : &DirEntry) -> bool {
    let filename = file.file_name().into_string().unwrap();
    if filename.ends_with("fna.gz") {
        return true;
    }
    else { 
        return false;
    }
}  // end of is_fasta_file




// opens parse fna files with needletail
// extracts records , filters out capsid and send 2 bits compressed sequence to some sketcher
fn process_file(file : &DirEntry)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let pathb = file.path();
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", file.file_name());
            std::process::exit(1);
        }
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        if strid.find("capsid").is_none() {
            // if we keep it we keep track of its id in file, we compress it with 2 bits per base
            let newseq = Sequence::new(&seqrec.seq(), 2);
            let seqwithid = IdSeq{rank : 0, id: strid, seq: newseq};
            to_sketch.push(seqwithid);
        }
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file



// TODO This function should have a version based on tokio::fs
// scan directory recursively, executing function cb.
// taken from fd_find
fn process_dir(dir: &Path, file_task: &dyn Fn(&DirEntry) -> Vec<IdSeq>, sender : &crossbeam_channel::Sender::<Vec<IdSeq>>) -> io::Result<SeqDict> {
    let mut nb_processed = 0;
    let mut seqdict = SeqDict::new(100000);
    // we checked that we have a directory
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            process_dir(&path, file_task, sender)?;
        } else {
            // check if entry is a fasta.gz file
            if is_fasta_file(&entry) {
                let mut to_sketch = file_task(&entry);
                // put a rank id in sequences;
                for i in 0..to_sketch.len() {
                    to_sketch[i].rank = nb_processed + i;
                    seqdict.0.push(to_sketch[i].id.clone());
                }
                nb_processed += to_sketch.len();
                // we must send to_sketch into channel to upper thread
                sender.send(to_sketch).unwrap();
            }
        }
    }
    println!("processed nb sequences : {}", nb_processed);
    //
    drop(sender);
    //
    // We must serialize and dump seqdict
    //
    Ok(seqdict)
}  // end of visit_dirs




// this function does the sketching and hnsw store of a whole directory
fn sketchandstore_dir(dirpath : &Path, sketcher_params : &SketcherParams, hnsw_params : &HnswParams) {
    // a queue of signature waiting to be inserted , size must be sufficient to benefit from threaded probminhash and insert
    let insertion_block_size = 1000;
    let mut insertion_queue : Vec<IdSeq>= Vec::with_capacity(insertion_block_size);
    // TODO must get ef_search from clap via hnswparams
    let hnsw = Hnsw::<u32, DistHamming>::new(hnsw_params.max_nb_conn , 700000, 16, hnsw_params.ef_search, DistHamming{});
    //
    // Sketcher allocation, we need reverse complement hashing
    //
    let kmer_revcomp_hash_fn = | kmer : &Kmer32bit | -> u32 {
        let canonical =  kmer.reverse_complement().min(*kmer);
        let hashval = probminhash::invhash::int32_hash(canonical.0);
        hashval
    };
    let sketcher = seqsketchjaccard::SeqSketcher::new(sketcher_params.kmer_size, sketcher_params.sketch_size);
    // to send IdSeq to sketch from reading thread to sketcher thread
    let (send, receive) = crossbeam_channel::bounded::<Vec<IdSeq>>(10_000);
    // launch process_dir in a thread or async
    crossbeam_utils::thread::scope(|scope| {
        let mut join_handles = Vec::new();
        // sequence sending, productor thread
        let sender_handle = scope.spawn(move |_|   {
            let _ = process_dir(dirpath, &process_file, &send);
            drop(send);
        });
        join_handles.push(sender_handle);
        // sequence reception, consumer thread
        let receptor_handle = scope.spawn(move |_| {
            // we must read messages, sketch and insert into hnsw
            let mut read_more = true;
            while read_more {
                // try read, if error is Disconnected we stop read and both threads are finished.
                let res_receive = receive.recv();
                match res_receive {
                    Err(RecvError) => { read_more = false;
                        // sketch the content of  insertion_queue if not empty
                        if insertion_queue.len() > 0 {
                            let sequencegroup_ref : Vec<&Sequence> = insertion_queue.iter().map(|s| &s.seq).collect();
                            let seq_id :  Vec<usize> = insertion_queue.iter().map(|s| s.rank).collect();
                            let signatures = sketcher.sketch_probminhash3a_kmer32bit(&sequencegroup_ref, kmer_revcomp_hash_fn);                            
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch for hnsw insertion
                            let mut data_for_hnsw = Vec::<(&Vec<u32>, usize)>::with_capacity(signatures.len());
                            for i in 0..signatures.len() {
                                data_for_hnsw.push((&signatures[i], seq_id[i]));
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
                            let sequencegroup_ref : Vec<&Sequence> = insertion_queue.iter().map(|s| &s.seq).collect();
                            let seq_id :  Vec<usize> = insertion_queue.iter().map(|s| s.rank).collect();
                            let signatures = sketcher.sketch_probminhash3a_kmer32bit(&sequencegroup_ref, kmer_revcomp_hash_fn);
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch, inserting unique id, for hnsw insertion
                            let mut data_for_hnsw = Vec::<(&Vec<u32>, usize)>::with_capacity(signatures.len());
                            for i in 0..signatures.len() {
                                data_for_hnsw.push((&signatures[i], seq_id[i]));
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
        }); // end of receptor thread
        join_handles.push(receptor_handle);
    }).unwrap();
    // We must dump hnsw to save "database"
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
            std::process::exit(1);
        }
        let dirpath = Path::new(&datadir);
        if !dirpath.is_dir() {
            println!("error not a directory : {:?}", datadir);
            std::process::exit(1);
        }
        // get sketching params
        let mut sketch_size = 8;
        if matches.is_present("size") {
            sketch_size = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("sketching size {}", sketch_size);
        }
        else {
            println!("using default sketch size {}", sketch_size);
        }
        //
        let mut kmer_size = 8;
        if matches.is_present("kmer_size") {
            kmer_size = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("kmer size {}", kmer_size);
        }
        else {
            println!("using default kmer size {}", kmer_size);
        }
        let sketch_params =  SketcherParams{kmer_size : kmer_size as usize, sketch_size : sketch_size as usize};  
        //
        let nbng;
        if matches.is_present("neighbours") {
            nbng = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("nb neighbours in hnsw size {}", nbng);
        }        
        else {
            std::process::exit(1);
        }
        // TODO pass parameters via clap
        let max_nb_conn = 48.min(3 * nbng as usize);
        let ef_search = 200;
        let hnswparams = HnswParams{nbng : nbng as usize, ef_search, max_nb_conn};
        //
        //
        sketchandstore_dir(&dirpath, &sketch_params, &hnswparams);


        //
 } // end of main