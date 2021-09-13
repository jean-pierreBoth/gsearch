//! tohnsw -d dir --ps size 
//! --ps gives the size of probminhash sketch

// must loop on sub directories , open gzipped files
// extracts complete genomes possiby many in one file (get rid of capsid records if any)
// compute probminhash sketch and store in a Hnsw.

// one thread should read sequences and do the probminhash
// another process should store in hnsw

// hnsw should also run in a query server mode after insertion.

 use clap::{App, Arg, SubCommand};


// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}


 fn main() {
    let _ = init_log();

 } // end of main