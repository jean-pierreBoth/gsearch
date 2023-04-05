
// GSEARCH v0.1.1
// Copyright 2021-2022, Jean Pierre Both and Jianshu Zhao.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.


//! Module gsearch
//!  
//! The commands admits 2 flags and 4 subcommands
//! 
//! The flags:
//! 
//!     --aa : set if data to process are Amino Acid sequences. Default is DNA.
//! 
//!     --pio nbfile : option to read compressed filesby blocks of n files then parallelize decompressing/fasta parsing. 
//!         Useful, with many cores if io lags behind hashing/hnsw insertion. to speed up io.  
//!         **Necessary to limit/custom the number of files or sequences simultanuously loaded in memory if files are very large (tens of Gb)**.
//! 
//! 
//! The 4 subcommands possible are to :
//! 
//!     - construct a hnsw database,  
//!  
//!     - adding elements in a already constructed hnsw database
//! 
//!     - search in hnsw database
//! 
//!     - ask neighbourhood statistics in the database and possibly ask for an embedding.
//!  
//! 1. ## subcommand  [--pio number] tohnsw --dir [-d] dir --sketch [-s] size --nbng [-n] nb --ef m [--seq]
//! 
//!     * \--dir : the name of directory containing tree of DNA files or Amino Acid files. 
//!   
//!     * \--sketch gives the size of probminhash sketch (integer value). Mandatory value.  
//! 
//!     * \--algo specifiy the sketching algorithm to be used. 
//! 
//!         SuperMinHash can specified by --algo super, --algo prob for asking ProbMinhash or --algo hll for hyperloglog sketching. 
//! 
//! 
//!     * \--kmer [-k] gives the size of kmer to use for generating probminhash (integer value). Mandatory argument. 
//!  
//!     * \--nbng [-n] gives the number of neihbours required in hnsw construction at each layer, in the range 24-64 is usual
//!             it doest not means you cannot ask for more neighbours in request.
//! 
//!     * \--ef optional integer value to optimize hnsw structure creation (default to 400)  
//! 
//!     * \--seq if we want a processing by sequences. Default is to concatenate all sequneces in a file
//!             in a large sequence.
//!  
//! 
//! 2. **sub command add : This command dedicated to adding new data to a hnsw structure.**.  
//!      The program reloads a previous dump of the hnsw structures.  
//!      tohnsw must (presently) be launched from the directory
//!      containing the dump as the program looks for the files "hnswdump.hnsw.data" and "hnswdump.hnsw.graph" created previously.  
//!      *In this case parameters corresponding to options --kmer  --sketch --nbng --ef and --algo are reloaded from file parameters.json*  .
//!      It is useless to pass them in command line.  
//!
//!
//! 3. **sub command request** 
//! 
//!  request --database [-b] basedirname --query [-r]  requestdir -n neighbours
//!     * \--database is the name of directory containing hnsw dump files and seqdict dump
//!     * \--requestdir is a directory containing list of fasta file containing sequence to search for
//!     * -n number of neighbours asked for. number of neighbours used in possible ann directive

//! 4. sub command ann

// must loop on sub directories , open gzipped files
// extracts complete genomes possiby many in one file (get rid of capsid records if any)
// compute probminhash sketch and store in a Hnsw.

// one thread should read sequences and do the probminhash
// another process should store in hnsw

// hnsw should also run in a query server mode after insertion.


//In fact as in basedirname there must be a file (processingparams.json) specifying sketch size and kmer size, these
// 2 options are useless in standard mode.

// We can use the same structure as the tohnsw module
// We parse a directory and send to a thread that do sketching and query
// we must enforce that the sketching size is the same as used in the database, so SketcherParams
// must have been dumped also in database directory.



use clap::{Arg, ArgMatches, Command, ArgAction};

use std::path::Path;

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

// our crate
use gsearch::dna::dnasketch::dna_process_tohnsw;
use gsearch::aa::aasketch::aa_process_tohnsw;
use gsearch::utils::*;

use kmerutils::sketcharg::{SketchAlgo, SeqSketcherParams};

//mod files;
use gsearch::dna::dnarequest;
use gsearch::aa::aarequest;

// parsing add command
// we need hnsw dir and directory containing new data, returns the 2-uple (hnsw_dir, newdata_dir)

#[derive(Clone, Debug)]
pub struct AddParams {
    pub(crate) hnsw_dir : String,
    //
    pub(crate) newdata_dir : String,
} // end of AddParams


impl AddParams {

    pub fn get_hnsw_dir(&self) -> &String { &self.hnsw_dir}

    pub fn get_newdata_dir(&self) -> &String { &self.newdata_dir}

}



//=======================================================================================


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
        log::error!("option --algo  super | prob | hll required");
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

    let data_t = if matches.contains_id("aa") {
        println!("aa option , processing of AA sequences");
        DataType::AA
    }
    else {
        println!("processing of DNA sequences");
        DataType::DNA
    };

    let sketch_params =  SeqSketcherParams::new(*kmer_size as usize, *sketch_size as usize, sketch_algo, data_t);  

    // max_nb_conn must be adapted to the number of neighbours we will want in searches.  
    // Maximum allowed nbng for hnswlib-rs is 256. Larger nbng will not work and default to 256.
    let max_nb_conn : u8 = 255.min(*nbng as u8);
    let hnswparams = HnswParams::new(1_500_000, *ef_construction, max_nb_conn);
    //
    let processing_params = ProcessingParams::new(hnswparams, sketch_params, block_processing);
    //
    log::debug!("subcommand parse_tohnsw : ok");
    //
    return Ok((datadir.clone(), processing_params));
} // end of parse_tohnsw


//===========================================================================================================



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


// Parsing request need hnsw database directory and data request dir. 
// All others args should be extracted from json reloads.
// The function returns (hnsw_dir, request_dir)

#[doc(hidden)]
fn parse_request_cmd(matches : &ArgMatches) -> Result<RequestParams, anyhow::Error> {
    log::debug!("in parse_request");
    //
    let to_check: &String = matches.get_one("request_dir").unwrap();  // as arg is required, we unwrap() !
    //
    let  request_dir = to_check.clone();
    let request_dirpath = Path::new(&to_check);
    if !request_dirpath.is_dir() {
        println!("error not a directory : {:?}", &request_dirpath);
        std::process::exit(1);
    }

    // parse database dir
    let to_check : &String = matches.get_one("database_dir").unwrap();
    let database_dir = to_check.clone();
    let database_dirpath = Path::new(&database_dir);
    if !database_dirpath.is_dir() {
        println!("error not a directory : {:?}", &database_dir);
        std::process::exit(1);
    }
    // decode nb_answers
    let nb_answers = *matches.get_one("nb_answers").unwrap();
    println!("nb answers {}", nb_answers);   
    //
    let req_params = RequestParams::new(database_dir.clone(), request_dir.clone(), nb_answers);
    return Ok(req_params);
}  // end of parse_request_cmd

//==========================================================================================

    // parse ann command

pub fn parse_ann_cmd( matches : &ArgMatches) -> Result<AnnParameters, anyhow::Error> {
    log::debug!("in parse_ann_cmd");
    let ask_stats = matches.get_flag("stats");
    let embed = matches.get_flag("embed");
    if ask_stats {
        log::info!("ann cmd with stats flag")
    }
    //
    if embed {
        log::info!("ann cmd with embed flag")
    }
    if !ask_stats && !embed {
        log::error!("ann cmd without flag");
        std::panic!("ann cmd without flag");
    }
    //
    let ann_params = AnnParameters::new(ask_stats, embed);
    //
    return Ok(ann_params);
} // end of parse_ann_cmd
 
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
        .arg(Arg::new("hnsw_dir")
            .short('d')
            .long("dir")
            .help("directory for storing database genomes")
            .required(true)
            .value_parser(clap::value_parser!(String))
            .action(ArgAction::Set)
            )
        .arg(Arg::new("kmer_size")
            .short('k')
            .long("kmer")
            .help("k-mer size to use")
            .required(true)
            .value_name("kmer_size")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
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
        .arg(Arg::new("neighbours")
            .short('n')
            .long("nbng")
            .help("Maximum allowed number of neighbors (M) in HNSW")
            .required(true)
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("ef_construct")
            .long("ef")
            .help("ef_construct in HNSW")
            .value_name("ef")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
        )
        .arg(Arg::new("sketch_algo")
            .required(true)
            .long("algo")
            .help("specifiy the algorithm to use for sketching: prob, super or hll")
            .value_parser(clap::value_parser!(String))
        )
        .arg(Arg::new("aa")                   // do we process amino acid file
            .long("aa")
            .value_name("AA")
            .help("Specificy amino acid processing")
        )
        .arg(Arg::new("seq")
            .long("seq")
            .help("--seq to get a processing by sequence")
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
    // ann command
    let ann_cmd = Command::new("ann")
        .about("annembed usage")
        .arg(Arg::new("stats")
            .long("stats")
            .short('s')
            .action(ArgAction::SetTrue)
            .help("to get stats on nb neighbours"))
        .arg(Arg::new("embed")
            .long("embed")
            .action(ArgAction::SetTrue)
            .help("--embed to do an embedding")
    );
    //
    // the global command
    //
    let matches = Command::new("gsearch")
        .version("0.1.1")
        .about("Approximate nearest neighbour search for microbial genomes based on minhash metric")
            .arg_required_else_help(true)
            .arg(Arg::new("pario")                // do we use parallel io
                .long("pio")
                .value_name("pio")
                .value_parser(clap::value_parser!(usize))
                .help("Parallel IO processing")
        )
        .subcommand(tohnsw_cmd)
        .subcommand(add_cmd)
        .subcommand(request_cmd)
        .subcommand(ann_cmd)
    .get_matches();

    // now we fill other parameters : parallel fasta parsing and adding mode in hnsw
    let nb_files_par: usize = *matches.get_one("pario").unwrap_or(&0usize);
    println!("parallel io, nb_files_par : {}", nb_files_par);
    //
    let hnsw_dir : String;
    let mut processing_params: Option<ProcessingParams> = None;
    //
    // treat tohnsw command
    //
    let cmd : CmdType;
    let mut add_params_opt : Option<AddParams> = None;
    let mut req_params_opt : Option<RequestParams> = None;
    let mut addseq = false;
    //
    // which command do we have
    //
    if let Some(tohnsw_match) = matches.subcommand_matches("tohnsw") {
        log::debug!("subcommand_matches got tohnsw command");
        let res = parse_tohnsw_cmd(tohnsw_match); 
        cmd = CmdType::TOHNSW;   
        match res {
            Ok(params) => { hnsw_dir = params.0;
                                                        processing_params = Some(params.1); },
            _          => { 
                            log::error!("parsing tohnsw command failed");
                            println!("exiting with error {}", res.as_ref().err().as_ref().unwrap());
                            log::error!("exiting with error {}", res.as_ref().err().unwrap());
                            std::process::exit(1);                                
                        },
        }
    } // end of tohnsw
    else if let Some(add_match) = matches.subcommand_matches("add") {
        // parsing add command
        log::debug!("subcommand_matches got add command"); 
        cmd = CmdType::ADD;
        addseq = true;
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
        hnsw_dir = add_params_opt.as_ref().unwrap().get_hnsw_dir().clone();
    }  // end of parsing add parameters
    else if let Some(req_match) = matches.subcommand_matches("request") {
        // parsing request parameters
        log::debug!("subcommand_matches got request command"); 
        cmd = CmdType::REQUEST;
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
        hnsw_dir = req_params_opt.as_ref().unwrap().get_hnsw_dir().clone();
    } // end of request
    else {
        log::error!("at least one command  tohnnsw, add or request must be given");
        std::panic!("at least one command  tohnnsw, add or request must be given");
    }

    let ann_params : AnnParameters;
    if let Some(ann_match) = matches.subcommand_matches("ann") {
        log::debug!("subcommand_matches got ann command"); 
        ann_params = parse_ann_cmd(ann_match).unwrap();
    }
    else {
        ann_params = AnnParameters::new(false, false)
    }
    //
    // Now we have parsed commands
    //
    let add_dir = if addseq { add_params_opt.as_ref().unwrap().get_hnsw_dir().clone() } else { String::from("")};
    let computing_params = ComputingParams::new(nb_files_par, addseq, add_dir);

    match cmd {
        CmdType::TOHNSW => {  // nothing to do we must have parsed before
                        }, 
             _          => {      // for case ADD or REQUEST we must reload
                           log::info!("reloading parameters from previous runs, in the current directory"); 
                            let cwd = std::env::current_dir().unwrap();
                            let reload_res = ProcessingParams::reload_json(&cwd);
                            if reload_res.is_ok()  {
                                processing_params = Some(reload_res.unwrap());
                            }
                            else {
                                std::panic!("cannot reload parameters (file parameters.json) from dir : {:?}", &cwd);
                            }
            },  // end of cae not TOHNSW
    };  // end of match

    let filter_params = FilterParams::new(0);
    if processing_params.is_none() {
        log::error!("error in processing command line, could not initilize ProcessingParams");
        std::panic!("error in processing command line, could not initilize ProcessingParams");
    }

    let processing_params = processing_params.unwrap();

    match cmd {
        CmdType::TOHNSW     => { treat_into_hnsw(&hnsw_dir, &filter_params, &processing_params, &computing_params); }

        CmdType::ADD        => { treat_into_hnsw(&hnsw_dir, &filter_params, &processing_params, &computing_params); }

        CmdType::REQUEST    => { 
                                    if req_params_opt.is_none() {
                                        log::error!("could not initialize request parameters");
                                        std::panic!("could not initialize request parameters");
                                    }
                                    treat_from_hnsw(req_params_opt.as_ref().unwrap(), &filter_params, 
                                                                &processing_params,
                                                                &computing_params,
                                                                &ann_params); }
    }
}  // end of main


//====================================================================

// function to create or add data into a hnsw database
fn treat_into_hnsw(hnswdir : &String, filter_params : &FilterParams, processing_params : &ProcessingParams, computing_parameters : &ComputingParams) {
     //
     let dirpath = Path::new(hnswdir);
     let data_type =  processing_params.get_sketching_params().get_data_t();
     match data_type {
         DataType::DNA => dna_process_tohnsw(&dirpath, &filter_params, &processing_params, &computing_parameters),
         DataType::AA => aa_process_tohnsw(&dirpath, &filter_params, &processing_params, &computing_parameters),
     }
} // end of treat_into_hnsw




// function to answer request from a hnsw database 
fn treat_from_hnsw(request_params : &RequestParams, filter_params : &FilterParams, processing_params : &ProcessingParams, 
        computing_params : &ComputingParams, ann_params : &AnnParameters) {

    let database_dir = request_params.get_hnsw_dir();
    let data_type =  processing_params.get_sketching_params().get_data_t();
    let ef_search = 5000;   // TODO ...
    // reload SeqDict
    let database_dirpath = Path::new(&database_dir);
    let seqname = database_dirpath.join("seqdict.json");
    log::info!("\n reloading sequence dictionary from {}", &seqname.display());
    let seqdict = SeqDict::reload_json(&seqname);
    let seqdict = match seqdict {
        Ok(seqdict ) => seqdict ,
        _ => {
            panic!("SeqDict reload from dump file  {} failed", seqname.display());
        }            
    };
    log::info!("reloading sequence dictionary from {} done", &seqname.display());
    // we have everything we want...
    match data_type {
        DataType::DNA => {
            if let Ok(mut seq_matcher) = dnarequest::get_sequence_matcher(request_params, &processing_params, 
                        &filter_params, &ann_params, &computing_params, &seqdict, ef_search) {
                if processing_params.get_block_flag() == false {
                    log::info!("sequence mode, trying to analyze..");
                    let _= seq_matcher.analyze();
                }
            }
            else {
                panic!("Error occurred in dnarequest::get_sequence_matcher");
            }
        }  // end DNA case
        DataType::AA => {
            if let Ok(mut seq_matcher) = aarequest::get_sequence_matcher(request_params, &processing_params, 
                        &filter_params, &ann_params, &computing_params, &seqdict,ef_search) {
                if processing_params.get_block_flag() == false {
                    log::info!("sequence mode, trying to analyze..");
                    let _= seq_matcher.analyze();
                }
            }
            else {
                panic!("Error occurred in aarequest::get_sequence_matcher");
            } 
        } // end of AA case
    }  // end of match on sequence type
    // 
} // end of treat_from_hnsw