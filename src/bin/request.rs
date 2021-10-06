//! Module request
//! try to match fasta sequence with repsect to database
//! 
//! request --database [-b] basedirname --query [-r]  requestdir --nbsearch [-n] nbanswers
//! 
//! - database is the name of directory containing hnsw dump files and seqdict dump
//! - requestdir is a file containing list of fasta file containing sequence to search for
//! 

// We can use the same structure as the tohnsw module
// We parse a directory and send to a thread that do sketching and query
// we must enforce that the sketching size is the same as used in the database, so SketcherParams
// must have been dumped also in database directory.
// TODO for now we pass it
// 



use clap::{App, Arg};

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

//
use std::time::{SystemTime};
use cpu_time::ProcessTime;


use std::path::{Path, PathBuf};
use std::fs::{OpenOptions, File};
use std::io::{Write,BufWriter};

use hnsw_rs::prelude::*;
use hnsw_rs::hnswio::{load_description, load_hnsw};
use kmerutils::base::{sequence::*, Kmer32bit, KmerT};
use kmerutils::sketching::*;
use kmerutils::sketching::seqsketchjaccard::*;

use archaea::utils::idsketch::{SeqDict, Id, IdSeq};
//mod files;
use archaea::utils::files::{process_dir,process_file};


// install a logger facility
pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}


/// An answer hs a rank (in which it is processed), the fasta Id corresponding to request sequence, and
/// the asked list of neighbours.
/// The neighbours are identified by an id in the database. To retrieve the fasta identity 
/// of the neighbour we will need the SeqDict
struct ReqAnswer<'a> {
    rank : usize,
    /// request id
    req_id : Id,
    ///
    neighbours : &'a Vec<Neighbour>,
}


impl <'a> ReqAnswer<'a> {
    pub fn new(rank : usize, req_id : Id, neighbours : &'a Vec<Neighbour>) -> Self {
        ReqAnswer { rank, req_id, neighbours}
    }

    /// dump answers to a File.
    fn dump(&self, seqdict : &SeqDict, out : &mut BufWriter<File>) -> std::io::Result<()> {
        // dump rank , fasta_id
        write!(out, "\n\n {} path {}, fasta_id {}", self.rank, self.req_id.get_path(), self.req_id.get_fasta_id())?;
        for n in self.neighbours {
            // get database identification of neighbour
            let database_id = seqdict.0[n.d_id].get_path();
            write!(out, "\n\t distance : {:.3E}  answer fasta id {}", n.distance, database_id)?;
        }
        Ok(())
    } // end of dump

}   // end of ReqAnswer



// TODO : possiblyy modify the way we pass requests via a whole directory.
/// this function does the sketching and hnsw store of a whole directory
/// dumpdir_path is the path to directory containing dump of Hnsw, database sequence dictionary and sketchparams.
/// request_dirpath is the directory containing the fasta files which are the requests.
/// 
fn sketch_and_request_dir(request_dirpath : &Path, hnsw : &Hnsw<u32,DistHamming>, seqdict : &SeqDict, 
                sketcher_params : &SeqSketcher, knbn : usize, ef_search : usize) {
    //
    log::trace!("sketch_and_request_dir processing dir {}", request_dirpath.to_str().unwrap());
    log::info!("sketch_and_request_dir {}", request_dirpath.to_str().unwrap());
    log::info!("sketch_and_request kmer size  {}  sketch size {} ", 
            sketcher_params.get_kmer_size(), sketcher_params.get_sketch_size());
    // creating an output file in the 
    let outname = "archea.answers";
    let outpath = PathBuf::from(outname.clone());
    let outfile = OpenOptions::new().write(true).create(true).truncate(true).open(&outpath);
    if outfile.is_err() {
        log::error!("SeqDict dump : dump could not open file {:?}", outpath.as_os_str());
        println!("SeqDict dump: could not open file {:?}", outpath.as_os_str());
        return Err("SeqDict Deserializer dump failed").unwrap();
    }    
    let mut outfile = BufWriter::new(outfile.unwrap());
    log::info!("dumping request answers in : {}", outname);
    //
    let start_t = SystemTime::now();
    let cpu_start = ProcessTime::now();
    // a queue of signature request , size must be sufficient to benefit from threaded probminhash and search
    let request_block_size = 500;
    let mut request_queue : Vec<IdSeq>= Vec::with_capacity(request_block_size);
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
            let res_nb_sent = process_dir(request_dirpath, &process_file, &send);
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
                            let sequencegroup_ref : Vec<&Sequence> = request_queue.iter().map(|s| s.get_sequence()).collect();
                            let seq_id :  Vec<Id> = request_queue.iter().map(|s| Id::new(s.get_path(), s.get_fasta_id())).collect();
                            let signatures = sketcher.sketch_probminhash3a_kmer32bit(&sequencegroup_ref, kmer_hash_fn);                            
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch for hnsw insertion
                            let knn_neighbours  = hnsw.parallel_search(&signatures, knbn, ef_search);
                            for i in 0..knn_neighbours.len() {
                                let answer = ReqAnswer::new(nb_request+i, seq_id[i].clone(), &knn_neighbours[i]);
                                if answer.dump(&seqdict, &mut outfile).is_err() {
                                    log::info!("could not dump answer for request id {}", answer.req_id.get_fasta_id());
                                }
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
                            let sequencegroup_ref : Vec<&Sequence> = request_queue.iter().map(|s| s.get_sequence()).collect();
                            // collect Id
                            let seq_id :  Vec<Id> = request_queue.iter().map(|s| Id::new(s.get_path(), s.get_fasta_id())).collect();
                            // computes hash signature
                            let signatures = sketcher.sketch_probminhash3a_kmer32bit(&sequencegroup_ref, kmer_hash_fn);
                            // we have Vec<u32> signatures we must go back to a vector of IdSketch, inserting unique id, for hnsw insertion
                            // parallel search
                            let knn_neighbours  = hnsw.parallel_search(&signatures, knbn, ef_search);
                            // construct and dump answers
                            for i in 0..knn_neighbours.len() {
                                let answer = ReqAnswer::new(nb_request+i, seq_id[i].clone(), &knn_neighbours[i]);
                                if answer.dump(&seqdict, &mut outfile).is_err() {
                                    log::info!("could not dump answer for request id {}", answer.req_id.get_fasta_id());
                                }
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
    }).unwrap();  // end of scope
    //
    let cpu_time = cpu_start.elapsed().as_secs();
    log::info!("process_dir : cpu time(s) {}", cpu_time);
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
    log::info!("process_dir : elapsed time(s) {}", elapsed_t);
} // end of sketch_and_request_dir 



/// reload hnsw from dump directory
/// We know filename : hnswdump.hnsw.data and hnswdump.hnsw.graph
fn reload_hnsw(dump_dirpath : &Path) -> Option<Hnsw<u32, DistHamming>> {
    // just concat dirpath to filenames and get pathbuf
    let graph_path = dump_dirpath.join("hnswdump.hnsw.graph");
    log::info!("reload_hnsw, loading graph from {}",graph_path.to_str().unwrap());
    let graphfile = OpenOptions::new().read(true).open(&graph_path);
    if graphfile.is_err() {
        println!("test_dump_reload: could not open file {:?}", graph_path.as_os_str());
        return None;
    }
    let mut graphfile = graphfile.unwrap();
    //
    let data_path = dump_dirpath.join("hnswdump.hnsw.data");
    log::info!("reload_hnsw, loading data from {}",data_path.to_str().unwrap());
    let datafile = OpenOptions::new().read(true).open(&data_path);
    if datafile.is_err() {
        println!("test_dump_reload: could not open file {:?}", data_path.as_os_str());
        return None;
    }
    let mut datafile = datafile.unwrap();
    //
    let start_t = SystemTime::now();
    let hnsw_description = load_description(&mut graphfile).unwrap();
    let hnsw : Hnsw<u32, DistHamming>= load_hnsw(&mut graphfile, &hnsw_description, &mut datafile).unwrap();
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
    log::info!("reload_hnsw : elapsed system time(s) {}", elapsed_t);
    //
    return Some(hnsw);
    //  
} // end of reload_hnsw





fn main() {
    let _ = init_log();
    //
    let matches = App::new("request")
        .arg(Arg::with_name("request_dir")
            .long("reqdir")
            .short("r")
            .takes_value(true)
            .help("name of directory containing request genomes to index"))
        .arg(Arg::with_name("database_dir")
            .long("datadir")
            .short("b")
            .takes_value(true)
            .help("name of directory containing database reference genomes"))            
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

        // decode matches, check for request_dir
        let request_dir;
        if matches.is_present("request_dir") {
            println!("decoding argument dir");
            request_dir = matches.value_of("request_dir").ok_or("").unwrap().parse::<String>().unwrap();
            if request_dir == "" {
                println!("parsing of request_dir failed");
                std::process::exit(1);
            }
        }
        else {
            println!("-r request_dir is mandatory");
            std::process::exit(1);
        }
        let request_dirpath = Path::new(&request_dir);
        if !request_dirpath.is_dir() {
            println!("error not a directory : {:?}", request_dirpath);
            std::process::exit(1);
        }

        // parse database dir
        let database_dir;
        if matches.is_present("database_dir") {
            println!("decoding argument dir");
            database_dir = matches.value_of("database_dir").ok_or("").unwrap().parse::<String>().unwrap();
            if database_dir == "" {
                println!("parsing of database_dir failed");
                std::process::exit(1);
            }
        }
        else {
            println!("-r database_dir is mandatory");
            std::process::exit(1);
        }
        let database_dirpath = Path::new(&database_dir);
        if !database_dirpath.is_dir() {
            println!("error not a directory : {:?}", database_dirpath);
            std::process::exit(1);
        }



        // get sketching params
        let mut sketch_size = 96;
        if matches.is_present("size") {
            sketch_size = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("do you know what you are doing, sketching size {}", sketch_size);
        }
        else {
            println!("will use dumped sketch size");
        }
        //
        let mut kmer_size = 8;
        if matches.is_present("kmer_size") {
            kmer_size = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("kmer size {}", kmer_size);
        }
        else {
            println!("will use dumped kmer size");
        }

        // in fact sketch_params must be initialized from the dump directory
        let _sketch_params =  SeqSketcher::new(kmer_size as usize, sketch_size as usize);  
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
        //
        // Do all dump reload, first sketch params. We reload smaller files first 
        // so that path errors are found early 
        //
        let sk_params = SeqSketcher::reload_json(database_dirpath);
        let sk_params = match sk_params {
            Ok(sk_params) => sk_params,
            _ => {
                panic!("SketchParams reload from dump dir {} failed", database_dirpath.to_str().unwrap());
            }
        };
        log::info!("sketch params reloaded kmer size : {}, sketch size {}", sk_params.get_kmer_size(), sk_params.get_sketch_size());
        // reload SeqDict
        let mut seqname = database_dir.clone();
        seqname.push_str("/seqdict.json");
        let seqdict = SeqDict::reload(&seqname);
        let seqdict = match seqdict {
            Ok(seqdict ) => seqdict ,
            _ => {
                panic!("SeqDict reload from dump file  {} failed", seqname);
            }            
        };
        // reload hnsw
        let hnsw = reload_hnsw(database_dirpath);
        let hnsw = match hnsw {
            Some(hnsw) => hnsw,
            _ => {
                panic!("hnsw reload from dump dir {} failed", database_dirpath.to_str().unwrap());
            }
        };

        // 
        let ef_search = 200;
        sketch_and_request_dir(&request_dirpath, &hnsw, &seqdict, &sk_params, nbng as usize, ef_search);

}  // end of main