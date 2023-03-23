// GSEARCH v0.1.1
// Copyright 2021-2022, Jean Pierre Both and Jianshu Zhao.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.


//! Module tohnsw 
//! gsearch tohnsw --dir [-d] dir --sketch [-s] size --nbng [-n] nb --ef m [--seq]
//! 
//! --dir : the name of directory containing tree of DNA files or Amino Acid files. 
//!   
//! --sketch gives the size of probminhash sketch (integer value). Mandatory value.  
//! 
//! --algo specifiy the sketching algorithm to be used. Default is ProbMinhash. SuperMinHash can specified by --algo super.
//!         passing --algo prob   for asking ProbMinhash is possible
//! 
//! --kmer [-k] gives the size of kmer to use for generating probminhash (integer value). Mandatory argument. 
//!  
//! --nbng [-n] gives the number of neihbours required in hnsw construction at each layer, in the range 24-64 is usual
//!             it doest not means you cannot ask for more neighbours in request.
//! 
//!  --ef optional integer value to optimize hnsw structure creation (default to 400)  
//! 
//!  --seq if we want a processing by sequences. Default is to concatenate all sequneces in a file
//!             in a large sequence.
//! 
//!  --aa : set if data to process are Amino Acid sequences. Default is DNA
//! 
//!  --pio : option to read compressed files and then parallelize decompressing/fasta parsing. 
//!         Useful, with many cores if io lags behind hashing/hnsw insertion. to speed up io.  
//!         **Necessary to limit/custom the number of files or sequences simultanuously loaded in memory if files are very large (tens of Gb)**.  
//! 
//! --add : This option is dedicated to adding new data to a hnsw structure.  
//!         The program reloads a previous dump of the hnsw structures. tohnsw must (presently) be launched from the directory
//!         containing the dump as the program looks for the files "hnswdump.hnsw.data" and "hnswdump.hnsw.graph" created previously.  
//!         **In this case parameters corresponding to options --kmer  --sketch --nbng --ef and --algo are reloaded from file parameters.json**.  
//!         It is useless to pass them in command line.

// must loop on sub directories , open gzipped files
// extracts complete genomes possiby many in one file (get rid of capsid records if any)
// compute probminhash sketch and store in a Hnsw.

// one thread should read sequences and do the probminhash
// another process should store in hnsw

// hnsw should also run in a query server mode after insertion.



//! Module request
//! try to match fasta sequence with repsect to database
//! 
//! gsearch request --database [-b] basedirname --query [-r]  requestdir -n neighbours \[ann\]
//! 
//! - database is the name of directory containing hnsw dump files and seqdict dump
//! - requestdir is a directory containing list of fasta file containing sequence to search for
//!  
//! -n number of neighbours asked for. number of neighbours used in possible ann directive
//! 
//! --pio : option to read compressed request files and then parallelize decompressing/fasta parsing. 
//!         Useful, with many cores if io lags behind hashing, to speed up io.  
//!         **The number of files or sequences simultanuously loaded in memory must be limited to fit in memory if files are very large (tens of Gb)**.  
//! 
//! 
//! --aa : set if data to process are Amino Acid sequences. Default is DNA
//!
//! \[ann\] is an optional subcommand asking for some statistics on distances between hnsw items or get an embedding of data
//!     --stats gives statistics on distances between neighbours
//!     --embed does a 2 dimensional embedding using the crate [annembed](https://crates.io/crates/annembed).
//! 
//! In fact as in basedirname there must be a file (processingparams.json) specifying sketch size and kmer size, these
//! 2 options are useless in standard mode.

// We can use the same structure as the tohnsw module
// We parse a directory and send to a thread that do sketching and query
// we must enforce that the sketching size is the same as used in the database, so SketcherParams
// must have been dumped also in database directory.



use clap::{Arg, ArgMatches, Command, arg};

use std::path::Path;

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

// our crate
use gsearch::dna::dnasketch::dna_process_tohnsw;
use gsearch::aa::aasketch::aa_process_tohnsw;
use gsearch::utils::*;

use kmerutils::sketcharg::{SketchAlgo, SeqSketcherParams};

//mod files;
use gsearch::utils::*;
use gsearch::dna::dnarequest;
use gsearch::aa::aarequest;


// install a logger facility
pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}


#[doc(hidden)]
fn parse_tohnsw(matches : &ArgMatches) -> Result<(), anyhow::Error> {
    log::debug!("in parse_tohnsw");
    //
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
    //
    // get sketching params
    //
    let sketch_size;
    if matches.is_present("sketch_size") {
        sketch_size = matches.value_of("sketch_size").ok_or("").unwrap().parse::<u16>().unwrap();
        println!("sketching size arg {}", sketch_size);
    }
    else {
        panic!("sketch_size is mandatory");
    }
    // sketching algorithm
    let mut sketch_algo = SketchAlgo::PROB3A;
    if matches.is_present("sketch_algo") {
        let algo = matches.value_of("sketch_algo").ok_or("").unwrap().parse::<String>().unwrap();
        println!("sketching algo {}", algo);
        if algo == String::from("super") {
            sketch_algo = SketchAlgo::SUPER;
        }
        else if algo == String::from("super2") {
            sketch_algo = SketchAlgo::SUPER2;
        }
        else if algo == String::from("prob") {
            sketch_algo = SketchAlgo::PROB3A;
        }
        else if algo == String::from("hll") {
            sketch_algo = SketchAlgo::HLL;
        }
        else {
            println!("unknown asketching algo");
            std::panic!("unknown sketching algo");
        }
    }
    // kmer size
    let kmer_size;
    if matches.is_present("kmer_size") {
        kmer_size = matches.value_of("kmer_size").ok_or("").unwrap().parse::<u16>().unwrap();
        println!("kmer size {}", kmer_size);
    }
    else {
        panic!(" kmer size is mandatory");
    }
    let sketch_params =  SeqSketcherParams::new(kmer_size as usize, sketch_size as usize, sketch_algo);  
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
    let mut ef_construction = 400;
    if matches.is_present("ef") {
        ef_construction = matches.value_of("ef").ok_or("").unwrap().parse::<usize>().unwrap();
        println!("ef construction parameters in hnsw construction {}", ef_construction);
    }
    else {
        println!("ef default used in construction {}", ef_construction);
    }
    // max_nb_conn must be adapted to the number of neighbours we will want in searches.  
    // Maximum allowed nbng for hnswlib-rs is 256. Larger nbng will not work and default to 256.
    let max_nb_conn : u8 = 255.min(nbng as u8);
    let hnswparams = HnswParams::new(1_500_000, ef_construction, max_nb_conn);
    //
    std::panic!("not yet");
} // end of parse_tohnsw


//============================================================================================



fn main() {
    let _ = init_log();

    let tohnsw_cmd  = Command::new("tohnsw")
        .about("Build HNSW graph database from a collection of database genomes based on minhash metric")
        .arg(Arg::with_name("directory")
            .short('d')
            .long("dir")
            .help("directory for storing database genomes")
            .required(true)
            .value_name("DIRECTORY")
            .takes_value(true),
            )
        .arg(Arg::with_name("kmer_size")
            .short('k')
            .long("kmer")
            .help("k-mer size to use")
            .required(true)
            .value_name("KMER_SIZE")
            .takes_value(true),
        )
        .arg(Arg::with_name("sketch_size")
            .short('s')
            .long("sketch")
            .help("sketch size of minhash to use")
            .required(true)
            .value_name("SKETCH_SIZE")
            .takes_value(true),
        )
        .arg(Arg::with_name("nbng")
            .short('n')
            .long("nbng")
            .help("Maximum allowed number of neighbors (M) in HNSW")
            .required(true)
            .value_name("NBNG")
            .takes_value(true),
        )
        .arg(Arg::with_name("ef_construct")
            .long("ef")
            .help("ef_construct in HNSW")
            .required(true)
            .value_name("EF")
            .takes_value(true),
        )
        .arg(Arg::with_name("sketch algo")
            .required(true)
            .help("specifiy the algorithm to use for sketching: prob, super or hll")
            .takes_value(true)
        )
        .arg(Arg::with_name("add")
            .long("add")
            .help("add to an existing HNSW index/database")
            .required(false)
            .value_name("ADD")
            .takes_value(false),
    );


    // add data to an already created hnsw. We need to specify only location of directory containing hnsw and directory containg data to add.
    // all others parameters are reloaded for json dumps.

    let add_cmd =  Command::new("add")
        .about("add file to a hnsw database")
        .arg(Arg::with_name("hnsw_dir")
            .required(true)
            .help("set the name of directory containing already constructed hnsw data")
            .takes_value(true)
        )
        .arg(Arg::with_name("add_dir")
            .required(true)
            .takes_value(true)
        );


    let request_cmd =  Command::new("request")
        .about("Request nearest neighbors of query genomes against a pre-built HNSW database/HNSW index")
        .arg(
            Arg::with_name("database_path")
            .short('b')
            .long("datadir")
            .value_name("DATADIR")
            .help("directory contains pre-built database files")
            .required(true)
            .takes_value(true),
        )
        .arg(
            Arg::with_name("nb_neighbors")
            .short('n')
            .long("nbanswer")
            .value_name("NEIGHBORS")
            .help("Sets the number of neighbors for the query")
            .required(true)
            .takes_value(true),
        )
        .arg(Arg::with_name("request_directory")
            .short('r')
            .long("request_directory")
            .value_name("REQUEST_DIRECTORY")
            .help("Sets the directory of request genomes")
            .required(true)
            .takes_value(true)
    );


    //
    // the global command
    //
    let matches = Command::new("gsearch")
        .version("0.1.1")
        .about("Approximate nearest neighbour search for microbial genomes based on minhash metric")
        .subcommand_required(true)
            .arg_required_else_help(true)
            .arg(Arg::with_name("aa")                   // do we process amino acid file
                .long("aa")
                .value_name("AA")
                .help("Specificy amino acid processing")
                .required(false)
                .takes_value(false),
            )
            .arg(Arg::with_name("pario")                // do we use parallel io
                    .long("pio")
                    .value_name("PIO")
                    .help("Parallel IO processing")
                    .required(false)
                    .takes_value(true),
            )
            .subcommand(tohnsw_cmd)
            .subcommand(request_cmd)
    .get_matches();

    // decode -aa and --pio options
    let data_type;
    if matches.is_present("aa") {
        println!("data to processs are AA data ");
        data_type = DataType::AA;
    }
    else {
        println!("data to processs are DNA data ");
        data_type = DataType::DNA;            
    }
    // now we fill other parameters : parallel fasta parsing and adding mode in hnsw
    let nb_files_par : usize;
    if matches.is_present("pario") {
        nb_files_par = matches.value_of("pario").ok_or("").unwrap().parse::<usize>().unwrap();
        println!("parallel io, nb_files_par : {}", nb_files_par);
    }
    else {
        nb_files_par = 0;
    }


/* 
    let matches = App::new("gsearch")
        .version("0.1.1")
        .author("Jean Pierre-Both <jeanpierre.both@gmail.com> and Jianshu Zhao <jianshu.zhao@gatech.edu>")
        .about("Approximate nearest neighbour search for microbial genomes based on minhash metric")
        .subcommand(
                    SubCommand::with_name("ann")
                        .about("An approximate nearest neighbor embedding with annembed")
                        .arg(
                            Arg::with_name("embed")
                                .long("embed")
                                .value_name("EMBEDDING")
                                .help("perform embedding")
                                .required(false),
                        )
                        .arg(
                            Arg::with_name("stats")
                                .long("stats")
                                .value_name("STATS")
                                .required(false)
                                .help("Displays stats for the ANN/HNSW graph")
                        ),
                ) // end ann subcommand
            )
        )
    ).get_matches();
    */
}