//! tohnsw -d dir --ps size 
//! --ps gives the size of probminhash sketch

// must loop on sub directories , open gzipped files
// extracts complete genomes possiby many in one file (get rid of capsid records if any)
// compute probminhash sketch and store in a Hnsw.

// one thread should read sequences and do the probminhash
// another process should store in hnsw

// hnsw should also run in a query server mode after insertion.

 use clap::{App, Arg, SubCommand};

 use std::io;
 use std::fs::{self, DirEntry};
 use std::path::Path;


// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
}


// returns true if file is a fasta file (possibly gzipped)
fn is_fasta_file(&DirEntry) -> bool {
    return false
}  // end of is_fasta_file



// scan directory recursively, executing function cb.
// taken from fd_find
fn visit_dirs(dir: &Path, cb: &dyn Fn(&DirEntry)) -> io::Result<()> {
    if dir.is_dir() {
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                visit_dirs(&path, cb)?;
            } else {
                // check if entry is a fasta.gz file
                cb(&entry);
            }
        }
    }
    Ok(())
}  // end of visit_dirs



 fn main() {
    let _ = init_log();

 } // end of main