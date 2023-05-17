
// GSEARCH v0.1.1
// Copyright 2021-2022, Jean Pierre Both and Jianshu Zhao.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.


//! # Module gsearch
//!  
//! ## Command and sub commands
//! 
//! The commands admits 2 flags and 4 subcommands
//! 
//! The flag:
//! 
//! * \--pio nbfile : option to read compressed filesby blocks of n files then parallelize decompressing/fasta parsing.   
//!         Useful, with many cores if io lags behind hashing/hnsw insertion. to speed up io.  
//!         **Necessary to limit/custom the number of files or sequences simultanuously loaded in memory if files are very large (tens of Gb)**.
//! 
//! The 4 subcommands possible are to :  
//! * construct a hnsw database,   
//! 
//! * add elements in a already constructed hnsw database  
//! 
//! * search in hnsw database  
//! 
//! * ask neighbourhood statistics in the database and possibly ask for an embedding.  
//!  
//! 1. ### subcommand  [--pio number] **tohnsw** [--aa]  --dir [-d] dir  --sketch [-s] size  --kmer kmersize  --algo name --nbng [-n] nb   --ef m [--seq]
//! 
//!     * \[--aa\] : set if data to process are Amino Acid sequences. Default is DNA.

//!     * \--dir  \[-d\]: the name of directory containing tree of DNA files or Amino Acid files. 
//!   
//!     * \--sketch gives the size of probminhash sketch (integer value). Mandatory value.  
//! 
//!     * \--algo specifiy the sketching algorithm to be used.   
//!         Use --algo super for superminhash , --algo prob for asking ProbMinhash or --algo hll for hyperloglog sketching.//! 
//! 
//!     * \--kmer [-k] gives the size of kmer to use for generating probminhash (integer value). Mandatory argument. 
//!  
//!     * \--nbng [-n] gives the number of neihbours required in hnsw construction at each layer, in the range 24-64 is usual
//!             it doest not means you cannot ask for more neighbours in request.
//! 
//!     * \--ef optional integer value to optimize hnsw structure creation (default to 400)  
//! 
//!     * \--seq if we want a processing by sequences. Default is to concatenate all sequences in a file
//!             in a large sequence.
//!  
//! 
//! 2. ### sub command **add**  --hnsw \[-b\] hnsw_dir --new \[-n\] directory
//! 
//!     * \--hnsw expects  the name of directory containing hnsw dump files and seqdict dump
//!     * \--new expects the name of the directory containing new data to add to the database.
//! 
//!      The command reloads a previous dump of the hnsw structures located in hnsw_dir and adds new data in it from directory
//!      containing the dump as the program looks for the files "hnswdump.hnsw.data" and "hnswdump.hnsw.graph" created previously.  
//!      *In this case parameters corresponding to options --aa --kmer  --sketch --nbng --ef and --algo are reloaded from file parameters.json*  .
//!      It is useless to pass them in command line.  
//!
//!
//! 3. ### sub command **request**  --hnsw \[-b\] databasedir --query [-r]  requestdir -n neighbours
//! 
//!     * \--hnsw expects the name of directory containing hnsw dump files and seqdict dump
//!     * \--query expects a directory containing list of fasta file containing sequence to search for
//!     * -n number of neighbours asked for, i.e number of neighbours asked for
//!
//! 4. ### sub command ann --hnsw \[\b\] hnsw_dir --embed or --stats (or both flags)
//! 
//!     * \--hnsw expects the name of directory containing hnsw dump files and seqdict dump, it is mandatory
//!     * \--embed    this flag ask for an embedding with default parameters. A csv file will be produced.
//!     * \--stats    this flag ask for neighbourhood statistics and node hubness histogram. See the crate [annembed](https://crates.io/crates/annembed)
//!     
//!     At least one of the two flags **embed** or **stats**  is required
//!     
//! 
//! ## Some hints on sketching to use
//! 
//!  The Probminhash algorithm takes multiplicity of kmers of kmers into account, so it is useful for cases where kmer multiplicity is important.
//!  It comes at a cost: if you have large files associated to large kmer the (hash) structure to store the count of billions of kmers needs much memory.  
//!  In this case consider using SuperMinHash or HyperLogLog which just record presence of kmers.  
//!  SuperMinHash in its present form sketch into f32 vectors but is quite fast. HyperLogLog is slower but sketch into u16 vector and can store billions
//!  of kmers. See the crate [probminhash](https://crates.io/crates/probminhash)
//!      
//! 
//! 



// TODO:
// - add a parsing of parameters to EmbedderParams in Ann subcommand
// - hnsw should also run in a query server mode after insertion.






use clap::{Arg, ArgMatches, Command, ArgAction};

use std::path::{Path, PathBuf};

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

    /// returns directory containing hnsw structiry containing present state of database
    pub fn get_hnsw_dir(&self) -> &String { &self.hnsw_dir}

    /// returns the directory containinf new data
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
    let ef_construction = matches.get_one("ef_construct").unwrap_or(&ef_construction_default);      
    println!("ef construction parameters in hnsw construction {}", ef_construction);
    
    let seq_flag = matches.get_flag("seq");
    if seq_flag {
        log::info!("setting seq flag true, block = false, will process files seq by seq");
    }
    else {
        log::info!("setting seq flag false, block = true, will process by concatenating seqs of file");
    }
    let block_processing = !seq_flag;
    //
    let data_t = if matches.contains_id("aa_opt") {
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
    if matches.contains_id("hnsw_dir") {
        println!("decoding argument dir");
        database_dir = matches.get_one("hnsw_dir").expect("");
        if database_dir == "" {
            println!("parsing of database_dir failed");
            std::process::exit(1);
        }
    }
    else {
        println!("-b database_dir is mandatory");
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
    let to_check : Option<&String> = matches.get_one("database_path");
    if to_check.is_none() {
        log::error!("no --query option passed");
        std::panic!("no --query option passed");
    }
    let to_check = to_check.unwrap();
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

pub fn parse_ann_cmd(matches : &ArgMatches) -> Result<AnnParameters, anyhow::Error> {
    log::debug!("in parse_ann_cmd");
    // parse database dir
    let database_dir : &String;
    if matches.contains_id("hnsw_dir") {
        println!("decoding argument dir");
        database_dir = matches.get_one("hnsw_dir").expect("");
        if database_dir == "" {
            println!("parsing of database_dir failed");
            std::process::exit(1);
        }
    }
    else {
        println!("-b database_dir is mandatory");
        std::process::exit(1);
    }
    //
    let ask_stats = matches.get_flag("stats");
    if ask_stats {
        log::info!("ann cmd with stats flag")
    }
    //
    let embed = matches.get_flag("embed");
    if embed {
        log::info!("ann cmd with embed flag")
    }
    // we must have at least one flag
    if !ask_stats && !embed {
        log::error!("ann cmd without flag");
        std::panic!("ann cmd without flag");
    }
    //
    let ann_params = AnnParameters::new(database_dir.clone(), ask_stats, embed);
    //
    return Ok(ann_params);
} // end of parse_ann_cmd
 
//============================================================================================

// we 
enum CmdType {
    TOHNSW,
    ADD,
    REQUEST,
    ANN,
}

fn main() {
    let _ = init_log();

    let tohnsw_cmd  = Command::new("tohnsw")
        .about("Build HNSW graph database from a collection of database genomes based on MinHash-like metric")
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
        .arg(Arg::new("aa_opt")                   // do we process amino acid file
            .long("aa")
            .num_args(0..=1)
            .default_missing_value("true")
            .help("Specificy amino acid processing, require no value")
        )
        .arg(Arg::new("seq")
            .long("seq")
            .action(clap::ArgAction::SetTrue)
            .help("--seq : sketching is done without concatenating sequences")
    );


    // add data to an already created hnsw. We need to specify only location of directory containing hnsw and directory containg data to add.
    // all others parameters are reloaded for json dumps.

    let add_cmd =  Command::new("add")
        .about("Add new genome files to a pre-built HNSW graph database")
        .arg(Arg::new("hnsw_dir")
            .required(true)
            .long("hnsw")
            .short('b')
            .value_parser(clap::value_parser!(String))
            .help("set the name of directory containing already constructed hnsw data")
        )
        .arg(Arg::new("newdata_dir")
            .required(true)
            .long("new")
            .short('n')
            .help("set directory containing new data")
            .value_parser(clap::value_parser!(String))
    );


    let request_cmd =  Command::new("request")
        .about("Request nearest neighbors of query genomes against a pre-built HNSW graph database/index")
        .arg(
            Arg::new("database_path")
            .short('b')
            .long("hnsw")
            .value_name("DATADIR")
            .help("directory contains pre-built database files")
            .required(true)
            .value_parser(clap::value_parser!(String))
        )
        .arg(
            Arg::new("nb_answers")
            .short('n')
            .long("nbanswers")
            .value_name("nb_answers")
            .help("Sets the number of neighbors for the query")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(usize))
            .required(true)
        )
        .arg(Arg::new("request_dir")
            .short('r')
            .long("query")
            .value_name("request_dir")
            .help("Sets the directory of request genomes")
            .value_parser(clap::value_parser!(String))
            .required(true)
    );
    // ann command
    let ann_cmd = Command::new("ann")
        .about("Approximate Nearest Neighbor Embedding using UMAP-like algorithm")
        .arg(Arg::new("hnsw_dir")
            .short('b')
            .long("hnsw")
            .help("directory containing hnsw")
            .required(true)
            .value_parser(clap::value_parser!(String))
            .action(ArgAction::Set)
        )
        .arg(Arg::new("stats")
            .long("stats")
            .short('s')
            .action(ArgAction::SetTrue)
            .help("to get stats on nb neighbours"))
        .arg(Arg::new("embed")
            .long("embed")
            .short('e')
            .action(ArgAction::SetTrue)
            .help("--embed to do an embedding")
    );
    //
    // the global command
    //
    let matches = Command::new("gsearch")
        .version("0.1.2")
        .about("Approximate nearest neighbour search for microbial genomes based on minhash metric")
            .arg_required_else_help(true)
            .arg(Arg::new("pario")                // do we use parallel io
                .long("pio")
                .value_name("pio")
                .value_parser(clap::value_parser!(usize))
                .help("Parallel IO processing"))
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
    let mut ann_params_opt : Option<AnnParameters> = None;
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
                                                        processing_params = Some(params.1);
                                                     },
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
    else if let Some(ann_match) = matches.subcommand_matches("ann") {
        // parsing ann command
        log::debug!("subcommand_matches got ann command");
        cmd = CmdType::ANN;
        let ann_params= parse_ann_cmd(ann_match).unwrap();
        hnsw_dir = ann_params.get_hnsw_dir().unwrap().clone();
        ann_params_opt = Some(ann_params);
    }
    else {
        log::error!("at least one command  tohnnsw, add, request or ann  must be given");
        std::panic!("at least one command  tohnnsw, add, request or ann must be given");
    }
    // 
    let _ann_params : AnnParameters;
    if ann_params_opt.is_none() {
        _ann_params = AnnParameters::default();
    }
    else{
        _ann_params = ann_params_opt.unwrap();
    }
    //
    // Now we have parsed commands
    //
    let add_dir = if addseq { add_params_opt.as_ref().unwrap().get_newdata_dir().clone() } else { String::from("")};
    if  addseq {
        log::info!("adding data from directory : {:?}", add_dir);
    }
    let computing_params = ComputingParams::new(nb_files_par, addseq, add_dir);

    match cmd {
        CmdType::TOHNSW => {  // nothing to do we must have parsed before
                        }, 
             _          => {      // for case ADD or REQUEST we must reload
                            log::info!("reloading parameters from previous runs, from directory : {:?}", &hnsw_dir); 
                            let hnsw_path = std::path::PathBuf::from(hnsw_dir.clone());
                            let reload_res = ProcessingParams::reload_json(&hnsw_path);
                            if reload_res.is_ok()  {
                                processing_params = Some(reload_res.unwrap());
                                log::info!("sketching parameters : {:?}", processing_params.as_ref().unwrap().get_sketching_params());
                                log::info!("block processing : {:?}", processing_params.as_ref().unwrap().get_block_flag());
                            }
                            else {
                                std::panic!("cannot reload parameters (file parameters.json) from dir : {:?}", &hnsw_path);
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
                                                                &computing_params); 
                                }  // end REQUEST

        CmdType::ANN            =>  {
                                    let hnsw_path = Path::new(&hnsw_dir);
                                    let type_name = reloadhnsw::get_hnsw_type(hnsw_path);
                                    if type_name.is_err() {
                                        log::error!("cannot get type of data in Hnsw");
                                        std::process::exit(1);
                                    }
                                    let type_name = type_name.unwrap();
                                    log::info!("got type in hnsw : {}", &type_name);
                                    #[cfg(any(feature="annembed_openblas-system", feature="annembed_openblas-static" , feature="annembed_intel-mkl"))]
                                    match type_name.as_str() {
                                        "u32" => {  // probminhash case
                                            let hnsw = reloadhnsw::reload_hnsw::<u32>(&hnsw_path).unwrap();
                                            if _ann_params.ask_stats() || _ann_params.embed() {
                                                log::info!("calling embed::get_graph_stats_embed");
                                                let _ = embed::get_graph_stats_embed(&hnsw, _ann_params.embed(), None);
                                            }
                                        }
                                        "u16" => {
                                            // hll case
                                            let hnsw = reloadhnsw::reload_hnsw::<u16>(&hnsw_path).unwrap();
                                            if _ann_params.ask_stats() || _ann_params.embed() {
                                                log::info!("calling embed::get_graph_stats_embed");
                                                let _ = embed::get_graph_stats_embed(&hnsw, _ann_params.embed(), None);
                                            }
                                        }
                                        "f32" => {
                                            // superminhash case
                                            let hnsw = reloadhnsw::reload_hnsw::<f32>(&hnsw_path).unwrap();
                                            if _ann_params.ask_stats() || _ann_params.embed() {
                                                log::info!("calling embed::get_graph_stats_embed");
                                                let _ = embed::get_graph_stats_embed(&hnsw, _ann_params.embed(), None);
                                            }
                                        }  
                                        _    => {
                                            log::error!("unknow type of data in hnsw, type : {}", &type_name);
                                        }                 
                                    } // end match
                                }  // end ANN 
    }
}  // end of main


//====================================================================

// function to create or add data into a hnsw database
fn treat_into_hnsw(hnswdir : &String, filter_params : &FilterParams, processing_params : &ProcessingParams, computing_parameters : &ComputingParams) {
     //
     let dirpath = PathBuf::from(hnswdir);
     let data_type =  processing_params.get_sketching_params().get_data_t();
     match data_type {
         DataType::DNA => dna_process_tohnsw(&dirpath, &filter_params, &processing_params, &computing_parameters),
         DataType::AA => aa_process_tohnsw(&dirpath, &filter_params, &processing_params, &computing_parameters),
     }
} // end of treat_into_hnsw




// function to answer request from a hnsw database 
fn treat_from_hnsw(request_params : &RequestParams, filter_params : &FilterParams, processing_params : &ProcessingParams, 
        computing_params : &ComputingParams) {

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
                        &filter_params, &computing_params, &seqdict, ef_search) {
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
                        &filter_params, &computing_params, &seqdict,ef_search) {
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
