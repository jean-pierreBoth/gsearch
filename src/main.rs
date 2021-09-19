//! tohnsw --dir [-d] dir --size [-s] size
//!  --dir : the name of directory containing tree of GCF and GCA files 
//! --size gives the size of probminhash sketch

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

// install a logger facility
fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");    
    return 1;
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
// extracts records , filters out capsid
fn process_file(file : &DirEntry) {
    let pathb = file.path();
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        // do we keep record ? we must get its id
        // if we keep it we keep track of its id in file, we sketch it
    }
} // end of process_file



// TODO This function should have a version bsed on tokio::fs
// scan directory recursively, executing function cb.
// taken from fd_find
fn process_dir(dir: &Path, cb: &dyn Fn(&DirEntry)) -> io::Result<()> {
    if dir.is_dir() {
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                process_dir(&path, cb)?;
            } else {
                // check if entry is a fasta.gz file
                if is_fasta_file(&entry) {
                    cb(&entry);
                }
            }
        }
    }
    Ok(())
}  // end of visit_dirs



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
            println!("dir argument is mandatory");
            std::process::exit(1);
        }
        //
        let mut sketch_size = 8;
        if matches.is_present("size") {
            sketch_size = matches.value_of("size").ok_or("").unwrap().parse::<u16>().unwrap();
            println!("sketching size {}", sketch_size);
        }
        else {
            println!("using default sketch size {}", sketch_size);
        }
        //
        let dirpath = Path::new(&datadir);
        //
        process_dir(dirpath, &process_file);
 } // end of main