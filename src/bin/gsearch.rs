
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



use clap::{Arg, ArgMatches, Command, ArgAction, arg};

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
fn parse_tohnsw_cmd(matches : &ArgMatches) -> Result<(String, ProcessingParams), anyhow::Error> {
    log::debug!("in parse_tohnsw");
    //
    let datadir : &String;
    if matches.contains_id("hnsw_dir") {
        println!("decoding argument dir");
        datadir = matches.get_one("hnsw_dir").expect("");
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
    let sketch_size: &usize;
    if matches.contains_id("sketch_size") {
        sketch_size = matches.get_one("sketch_size").expect("required");
        println!("sketching size arg {}", sketch_size);
    }
    else {
        panic!("sketch_size is mandatory");
    }
    // sketching algorithm
    let sketch_algo : SketchAlgo;
    if matches.contains_id("sketch_algo") {
        let algo_opt : Option<&String> = matches.get_one("sketch_algo");
        if algo_opt.is_none() {
            std::panic!("algo needs a choice");  
        }
        let algoname = algo_opt.unwrap().clone();
        println!("sketching algo {}", algoname);
        if algoname == String::from("super") {
            sketch_algo = SketchAlgo::SUPER;
        }
        else if algoname == String::from("super2") {
            sketch_algo = SketchAlgo::SUPER2;
        }
        else if algoname == String::from("prob") {
            sketch_algo = SketchAlgo::PROB3A;
        }
        else if algoname == String::from("hll") {
            sketch_algo = SketchAlgo::HLL;
        }
        else {
            println!("unknown asketching algo");
            std::panic!("unknown sketching algo");
        }
    }
    else {
        std::panic!("sketching algo required");
    }
    // kmer size
    let kmer_size : &usize;
    if matches.contains_id("kmer_size") {
        kmer_size = matches.get_one("kmer_size").expect("required");
        println!("kmer size {}", kmer_size);
    }
    else {
        panic!(" kmer size is mandatory");
    }
    let sketch_params =  SeqSketcherParams::new(*kmer_size as usize, *sketch_size as usize, sketch_algo);  
    //
    let nbng: &usize;
    if matches.contains_id("neighbours") {
        nbng = matches.get_one("neighbours").expect("required");
        println!("nb neighbours you will need in hnsw requests {}", nbng);
    }        
    else {
        println!("-n nbng is mandatory");
        std::process::exit(1);
    }
    //
    let ef_construction_default = 400usize;
    let ef_construction = matches.get_one("ef").unwrap_or(&ef_construction_default);      
    println!("ef construction parameters in hnsw construction {}", ef_construction);
    
    let block_processing = if matches.contains_id("seq") {
        println!("seq option , will process by concatenating and splitting ");
        false
    }
    else {
        println!("seq option , will process by concatenating");
        true
    };
    // max_nb_conn must be adapted to the number of neighbours we will want in searches.  
    // Maximum allowed nbng for hnswlib-rs is 256. Larger nbng will not work and default to 256.
    let max_nb_conn : u8 = 255.min(*nbng as u8);
    let hnswparams = HnswParams::new(1_500_000, *ef_construction, max_nb_conn);
    //
    let processing_params = ProcessingParams::new(hnswparams, sketch_params, block_processing);
    return Ok((datadir.clone(), processing_params));
} // end of parse_tohnsw


//===========================================================================================================

// parsing add command
// we need hnsw dir and directory containing new data, returns the 2-uple (hnsw_dir, newdata_dir)

#[derive(Clone, Debug)]
struct AddParams {
    pub(crate) hnsw_dir : String,
    //
    pub(crate) newdata_dir : String,
} // end of AddParams


#[doc(hidden)]
fn parse_add_cmd(matches : &ArgMatches) -> Result<AddParams, anyhow::Error> {
    log::debug!("in add_tohnsw");
    //
    // parse database dir
    let database_dir : &String;
    if matches.contains_id("database_dir") {
        println!("decoding argument dir");
        database_dir = matches.get_one("hnsw_dir").expect("");
        if database_dir == "" {
            println!("parsing of database_dir failed");
            std::process::exit(1);
        }
    }
    else {
        println!("-r database_dir is mandatory");
        std::process::exit(1);
    }
    //
    // parse new data directory
    //
    let newdata_dir: &String;
    if matches.contains_id("newdata_dir") {
        println!("decoding argument dir");
        newdata_dir = matches.get_one("newdata_dir").expect("");
        if newdata_dir == "" {
            println!("parsing of newdata_dir failed");
            std::process::exit(1);
        }
    }
    else {
        println!("-a newdata_dir is mandatory");
        std::process::exit(1);
    }    
    //
    let add_params =  AddParams { hnsw_dir : database_dir.clone(), newdata_dir : newdata_dir.clone()};
    //
    return Ok(add_params);
} // end of parse_add_cmd

//============================================================================================================

struct RequestParams {
    pub(crate) hnsw_dir : String,
    //
    pub(crate) req_dir : String,
    //
    pub(crate) nb_answers : usize,
} // end of RequestParams


// Parsing request need hnsw database directory and data request dir. 
// All others args should be extracted from json reloads.
// The function returns (hnsw_dir, request_dir)

#[doc(hidden)]
fn parse_request_cmd(matches : &ArgMatches) -> Result<RequestParams, anyhow::Error> {
    log::debug!("in parse_request");
    //
    let request_dir : &String;
    if matches.contains_id("request_dir") {
        println!("decoding argument dir");
        request_dir = matches.get_one("request_dir").unwrap();  // as arg is required, we unwrap() !
    }
    //
    let to_check = request_dir.clone();
    let request_dirpath = Path::new(&to_check);
    if !request_dirpath.is_dir() {
        println!("error not a directory : {:?}", &request_dirpath);
        std::process::exit(1);
    }

    // parse database dir
    let database_dir : &String;
    if matches.contains_id("database_dir") {
        println!("decoding argument dir");
        database_dir = matches.get_one("database_dir").unwrap();
    }
    let to_check = database_dir.clone();
    let database_dirpath = Path::new(&to_check);
    if !database_dirpath.is_dir() {
        println!("error not a directory : {:?}", &to_check);
        std::process::exit(1);
    }
    // decode nb_answers
    let nb_answers : usize;
    if matches.contains_id("nb_answers") {
        nb_answers = *matches.get_one("nb_answers").unwrap();
        println!("nb answers {}", nb_answers);
    }    
    //
    let req_params = RequestParams{hnsw_dir: database_dir.clone(), req_dir :  request_dir.clone(), nb_answers};
    return Ok(req_params);
}  // end of parse_request_cmd



//============================================================================================

// we 
enum CmdType {
    TOHNSW,
    ADD,
    REQUEST,
}

fn main() {
    let _ = init_log();

    let tohnsw_cmd  = Command::new("tohnsw")
        .about("Build HNSW graph database from a collection of database genomes based on minhash metric")
        .arg(Arg::new("hnswdir")
            .short('d')
            .long("dir")
            .help("directory for storing database genomes")
            .required(true)
            .value_parser(clap::value_parser!(String))
            .action(ArgAction::Set)
            .value_name("hnswdir")
            )
        .arg(Arg::new("kmer_size")
            .short('k')
            .long("kmer")
            .help("k-mer size to use")
            .required(true)
            .value_name("kmer_size")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(u16))
        )
        .arg(Arg::new("sketch_size")
            .short('s')
            .long("sketch")
            .help("sketch size of minhash to use")
            .required(true)
            .value_name("sketch_size")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("nbng")
            .short('n')
            .long("nbng")
            .help("Maximum allowed number of neighbors (M) in HNSW")
            .required(true)
            .value_name("nbng")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("ef_construct")
            .long("ef")
            .help("ef_construct in HNSW")
            .required(true)
            .value_name("ef")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("sketch algo")
            .required(true)
            .help("specifiy the algorithm to use for sketching: prob, super or hll")
            .value_parser(clap::value_parser!(String))
        )
        .arg(Arg::new("seq")
            .long("seq")
            .help("--seq to get a processing by sequence"))
        .arg(Arg::new("add")
            .long("add")
            .help("add to an existing HNSW index/database")
            .required(false)
            .value_name("add")
    );


    // add data to an already created hnsw. We need to specify only location of directory containing hnsw and directory containg data to add.
    // all others parameters are reloaded for json dumps.

    let add_cmd =  Command::new("add")
        .about("add file to a hnsw database")
        .arg(Arg::new("hnsw_dir")
            .required(true)
            .value_parser(clap::value_parser!(String))
            .help("set the name of directory containing already constructed hnsw data")
        )
        .arg(Arg::new("add_dir")
            .required(true)
            .help("set directory containing new data")
            .value_parser(clap::value_parser!(String))
        );


    let request_cmd =  Command::new("request")
        .about("Request nearest neighbors of query genomes against a pre-built HNSW database/HNSW index")
        .arg(
            Arg::new("database_path")
            .short('b')
            .long("datadir")
            .value_name("DATADIR")
            .help("directory contains pre-built database files")
            .required(true)
            .value_parser(clap::value_parser!(String))
        )
        .arg(
            Arg::new("nb_answers")
            .short('n')
            .long("nbanswer")
            .value_name("nb_answers")
            .help("Sets the number of neighbors for the query")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
            .required(true)
        )
        .arg(Arg::new("request_dir")
            .short('r')
            .long("request_directory")
            .value_name("request_dir")
            .help("Sets the directory of request genomes")
            .value_parser(clap::value_parser!(String))
            .required(true)
    );


    //
    // the global command
    //
    let matches = Command::new("gsearch")
        .version("0.1.1")
        .about("Approximate nearest neighbour search for microbial genomes based on minhash metric")
            .arg_required_else_help(true)
            .arg(Arg::new("aa")                   // do we process amino acid file
                .long("aa")
                .value_name("AA")
                .help("Specificy amino acid processing")
                .required(false)
            )
            .arg(Arg::new("pario")                // do we use parallel io
                .long("pio")
                .value_name("pio")
                .value_parser(clap::value_parser!(usize))
                .help("Parallel IO processing")
        )
        .subcommand(tohnsw_cmd)
        .subcommand(add_cmd)
        .subcommand(request_cmd)
    .get_matches();

    // decode -aa and --pio options
    let data_type;
    if matches.contains_id("aa") {
        println!("data to processs are AA data ");
        data_type = DataType::AA;
    }
    else {
        println!("data to processs are DNA data ");
        data_type = DataType::DNA;            
    }
    // now we fill other parameters : parallel fasta parsing and adding mode in hnsw
    let nb_files_par : usize;
    if matches.contains_id("pario") {
        let nb_files_par: usize = *matches.get_one("pario").unwrap_or(&0usize);
        println!("parallel io, nb_files_par : {}", nb_files_par);
    }
    else {
        nb_files_par = 0;
    }
    //
    let hnsw_dir : String;
    let processing_params;
    let hnsw_params : HnswParams;
    let seqsketch_params : SeqSketcherParams;
    //
    // treat tohnsw command
    //
    let mut nb_cmd = 0;
    let cmd : CmdType;

    if let Some(tohnsw_match) = matches.subcommand_matches("tohnsw") {
        nb_cmd += 1;
        log::debug!("subcommand_matches got tohnsw command");
        let res = parse_tohnsw_cmd(tohnsw_match); 
        cmd = CmdType::TOHNSW;   
        match res {
            Ok(params) => { hnsw_dir = params.0;
                                                        processing_params = params.1; },
            _          => { 
                            log::error!("parsing tohnsw command failed");
                            println!("exiting with error {}", res.as_ref().err().as_ref().unwrap());
                            log::error!("exiting with error {}", res.as_ref().err().unwrap());
                            std::process::exit(1);                                
                        },
        }
    } // end of tohnsw

    let mut addseq = false;
    let add_params_opt : Option<AddParams> = None;
    if let Some(add_match) = matches.subcommand_matches("add") {
        log::debug!("subcommand_matches got add command"); 
        if nb_cmd >= 1 {
            log::error!("only one main subcommand is possible, tohnsw, add or request");
            std::panic!("only one main subcommand is possible, tohnsw, add or request");
        }
        else {
            cmd = CmdType::ADD;
            addseq = true;
        }
        nb_cmd += 1;
        let res = parse_add_cmd(add_match);
        match res {
            Ok(params) => {
                            add_params_opt = Some(params);
                        }
            _          => {
                            log::error!("parsing add command failed");
                            println!("exiting with error {}", res.as_ref().err().as_ref().unwrap());
                            log::error!("exiting with error {}", res.as_ref().err().unwrap());
                            std::process::exit(1);                  
                        }
        }  
    }  // end of parsing add parameters

    // parsing request parameters

    let req_params_opt : Option<RequestParams> = None;
    if let Some(req_match) = matches.subcommand_matches("request") {
        log::debug!("subcommand_matches got request command"); 
        if nb_cmd >= 1 {
            log::error!("only one main subcommand is possible, tohnsw, add or request");
            std::panic!("only one main subcommand is possible, tohnsw, add or request");
        }
        else {
            cmd = CmdType::REQUEST;
        }
        nb_cmd += 1;
        let res = parse_request_cmd(req_match);
        match res {
            Ok(params) => { 
                            req_params_opt = Some(params);
            }
           _           => {
                            log::error!("parsing request command failed");
                            println!("exiting with error {}", res.as_ref().err().as_ref().unwrap());
                            log::error!("exiting with error {}", res.as_ref().err().unwrap());
                            std::process::exit(1);                  
                        }
        }
    } // end of request

    //
    // Now we have parsed commands
    //
    let computing_params = ComputingParams::new(nb_files_par, addseq);


    match cmd {
        CmdType::TOHNSW => {  // nothing to do we must have parsed before
                        }, 
             _          => {      // for case ADD or REQUEST we must reload
                           log::info!("reloading parameters from previous runs, in the current directory"); 
                            let cwd = std::env::current_dir().unwrap();
                            let reload_res = ProcessingParams::reload_json(&cwd);
                                if reload_res.is_ok()  {
                                    processing_params = reload_res.unwrap()
                                }
                            else {
                                std::panic!("cannot reload parameters (file parameters.json) from dir : {:?}", &cwd);
                            }
            },  // end of cae not TOHNSW
    };  // end of match


    match cmd {
        CmdType::TOHNSW     => { treat_into_hnsw(); }
        CmdType::ADD        => { treat_into_hnsw(); }
        CmdType::REQUEST    => { treat_from_hnsw(); }
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
                            Arg::new("embed")
                                .long("embed")
                                .value_name("EMBEDDING")
                                .help("perform embedding")
                                .required(false),
                        )
                        .arg(
                            Arg::new("stats")
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
}  // end of main


//====================================================================

// function to create or add data into a hnsw database
fn treat_into_hnsw(hnswdir : &String, filter_params : &FilterParams, processing_parameters : &ProcessingParams, other_params : ComputingParams) {
    // do not filter small seqs when running file in a whole block
     let filter_params = FilterParams::new(0);
     //
     match data_type {
         DataType::DNA => dna_process_tohnsw(&dirpath, &filter_params, &processing_parameters, &other_params),
         DataType::AA => aa_process_tohnsw(&dirpath, &filter_params, &processing_parameters, &other_params),
     }
} // end of treat_into_hnsw


// function to answer request from a hnsw database 
fn treat_from_hnsw() {}