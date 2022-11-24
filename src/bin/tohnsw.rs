// ARCHAEA v0.1.0
// Copyright 2021-2022, Jean Pierre Both and Jianshu Zhao.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.




//! tohnsw --dir [-d] dir --sketch [-s] size --nbng [-n] nb --ef m [--seq]
//! 
//! --dir : the name of directory containing tree of DNA files or Amino Acid files. 
//!   
//! --sketch gives the size of probminhash sketch (integer value). Mandatory value.  
//! 
//! --kmer [-k] gives the size of kmer to use for generating probminhash (integer value). Mandatory argument. 
//!  
//! --nbng [-n] gives the number of neihbours required in hnsw construction at each layer, in the range 24-64 is usual
//!             it doest not means you cannot ask for more neighbours in request.
//! 
//!  -- ef optional integer value to optimize hnsw structure creation (default to 400)  
//! 
//!  --seq if we want a processing by sequences. Default is to concatenate all sequneces in a file
//!             in a large sequence.
//! 
//!  --aa : set if data to process are Amino Acid sequences. Default is DNA
//! 
//!  --pio : option to read uncompressed files and then parallelize decompressing/fasta parsing. 
//!          Useful, with many cores if io lags behind hashing/hnsw insertion.
//!          Implemented only for dna files presently

// must loop on sub directories , open gzipped files
// extracts complete genomes possiby many in one file (get rid of capsid records if any)
// compute probminhash sketch and store in a Hnsw.

// one thread should read sequences and do the probminhash
// another process should store in hnsw

// hnsw should also run in a query server mode after insertion.

use clap::{Command, Arg};

use std::path::Path;

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

//


// our crate

use gsearch::dna::dnasketch::dna_process_tohnsw;
use gsearch::aa::aasketch::aa_process_tohnsw;
use gsearch::utils::*;

//=========================================================================

// install a logger facility
pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}

// this function does the sketching in small kmers (less than 14 bases) and hnsw store of a whole directory, version before generic one

fn main() {
    let _ = init_log();
    //
    //
    let matches = Command::new("tohnsw")
        .arg(Arg::new("dir")
            .long("dir")
            .short('d')
            .takes_value(true)
            .required(true)
            .help("name of directory containing genomes to index"))
        .arg(Arg::new("kmer_size")
            .long("kmer")
            .short('k')
            .required(true)
            .takes_value(true)
            .help("expecting a kmer size"))
        .arg(Arg::new("sketch_size")
            .long("sketch")
            .short('s')
            .required(true)
            .takes_value(true)
            .help("size of probinhash sketch, default to 256"))
        .arg(Arg::new("neighbours")
            .long("nbng")
            .short('n')
            .required(true)
            .takes_value(true)
            .help("must specify number of neighbours in hnsw"))
        .arg(Arg::new("ef")
            .long("ef")
            .default_value("400")
            .help("parameters neighbour search at creation"))
        .arg(Arg::new("aa")
            .help("to specify amino acid seq processing")
            .long("aa")
            .takes_value(false))
        .arg(Arg::new("seq")
            .long("seq")
            .takes_value(false)
            .help("--seq to get a processing by sequence"))
        .arg(Arg::new("pario")
            .long("pio")
            .takes_value(false)
            .help("--pio to optimize io"))
        .get_matches();
    //
    // by default we process DNA files in one large sequence block
    let mut block_processing = true;
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
        let sketch_size;
        if matches.is_present("sketch_size") {
            sketch_size = matches.value_of("sketch_size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("sketching size arg {}", sketch_size);
        }
        else {
            panic!("sketch_size is mandatory");
        }
        //
        let kmer_size;
        if matches.is_present("kmer_size") {
            kmer_size = matches.value_of("kmer_size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("kmer size {}", kmer_size);
        }
        else {
            panic!(" kmer size is mandatory");
        }
        let sketch_params =  SeqSketcherParams::new(kmer_size as usize, sketch_size as usize);  
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
        // do we use block processing, recall that default is yes
        if matches.is_present("seq") {
            println!("seq option , will process every sequence independantly ");
            block_processing = false;
        }
        //
        let data_type;
        if matches.is_present("aa") {
            println!("data to processs are AA data ");
            data_type = DataType::AA;
        }
        else {
            println!("data to processs are DNA data ");
            data_type = DataType::DNA;            
        }
        // parallel io
        let pario: bool;
        if matches.is_present("pario") {
            println!("parallel io");
            pario = true;
        }
        else {
            pario = false;
        }
        // We have everything   
        // max_nb_conn must be adapted to the number of neighbours we will want in searches.
        
        // Maximum allowed nbng for hnswlib-rs is 256. Larger nbng will not work and default to 256.
        let max_nb_conn : u8 = 255.min(nbng as u8);
        let hnswparams = HnswParams::new(1_500_000, ef_construction, max_nb_conn);
        //
        // do not filter small seqs when running file in a whole block
        let filter_params = FilterParams::new(0);
        let processing_parameters = ProcessingParams::new(hnswparams, sketch_params, block_processing);
        //
        match data_type {
            DataType::DNA => dna_process_tohnsw(&dirpath, &filter_params, &processing_parameters, pario),
            DataType::AA => aa_process_tohnsw(&dirpath, &filter_params, &processing_parameters),
        }
        //
 } // end of main