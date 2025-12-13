use std::fs::File;
use std::io::{self, BufRead};

use clap::{Arg, Command};
use env_logger::Builder;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use sourmash::_hash_murmur;
use sourmash::encodings::HashFunctions;
use sourmash::sketch::minhash::KmerMinHash;

pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    1
}

fn main() -> io::Result<()> {
    let _ = init_log();

    let matches = Command::new("superaai")
        .version("0.3.2")
        .about("Compute Average Amino Acid Identity (AAI) via FracMinHash/Sourmash for genomes")
        .arg(
            Arg::new("query_list")
                .short('q')
                .long("ql")
                .value_name("FILE")
                .help("File containing list of query protein paths (.faa format, .gz supported)")
                .required(true),
        )
        .arg(
            Arg::new("ref_list")
                .short('r')
                .long("rl")
                .value_name("FILE")
                .help(
                    "File containing list of reference protein paths (.faa format, .gz supported)",
                )
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output file to write results")
                .required(true),
        )
        .arg(
            Arg::new("kmer_size")
                .short('k')
                .long("kmer")
                .value_name("INT")
                .help("K-mer size for MinHash calculation")
                .required(false)
                .default_value("7"),
        )
        .arg(
            Arg::new("scaled")
                .short('l')
                .long("scaled")
                .value_name("INT")
                .help("Scaled factor for MinHash calculation")
                .required(false)
                .default_value("100"),
        )
        .arg(
            Arg::new("sketch")
                .short('s')
                .long("sketch")
                .value_name("INT")
                .help("Sketch size for MinHash (number of hashes)")
                .required(false)
                .default_value("5120"),
        )
        .get_matches();

    let query_list_path = matches.get_one::<String>("query_list").unwrap();
    let ref_list_path = matches.get_one::<String>("ref_list").unwrap();
    let output_file_path = matches.get_one::<String>("output").unwrap();
    let kmer_size: u32 = matches
        .get_one::<String>("kmer_size")
        .unwrap()
        .parse()
        .expect("Invalid k-mer size");
    let scaled: u32 = matches
        .get_one::<String>("scaled")
        .unwrap()
        .parse()
        .expect("Invalid scaled factor");
    let sketch: u32 = matches
        .get_one::<String>("sketch")
        .unwrap()
        .parse()
        .expect("Invalid sketch size");

    let query_files: Vec<String> = {
        let file = File::open(query_list_path)?;
        io::BufReader::new(file)
            .lines()
            .filter_map(io::Result::ok)
            .collect()
    };

    let ref_files: Vec<String> = {
        let file = File::open(ref_list_path)?;
        io::BufReader::new(file)
            .lines()
            .filter_map(io::Result::ok)
            .collect()
    };

    let results: Vec<String> = query_files
        .par_iter()
        .flat_map(|query_file| {
            ref_files.par_iter().map(move |ref_file| {
                // Computing MinHash for query file
                let mut query_mh = KmerMinHash::new(
                    scaled,
                    kmer_size,
                    HashFunctions::Murmur64Protein,
                    42,
                    false,
                    sketch,
                );
                let mut query_reader =
                    parse_fastx_file(query_file).expect("Cannot open query file");
                while let Some(record) = query_reader.next() {
                    let record = record.expect("Failed to read query record");
                    // Process each k-mer from the sequence
                    for kmer in record.seq().windows(kmer_size as usize) {
                        let hash = _hash_murmur(kmer, 42); // Adjust the hash function as needed
                        query_mh.add_hash(hash);
                    }
                }

                // Computing MinHash for reference file
                let mut ref_mh = KmerMinHash::new(
                    scaled,
                    kmer_size,
                    HashFunctions::Murmur64Protein,
                    42,
                    false,
                    sketch,
                );
                let mut ref_reader =
                    parse_fastx_file(ref_file).expect("Cannot open reference file");
                while let Some(record) = ref_reader.next() {
                    let record = record.expect("Failed to read reference record");
                    // Process each k-mer from the sequence
                    for kmer in record.seq().windows(kmer_size as usize) {
                        let hash = _hash_murmur(kmer, 42); // Adjust the hash function as needed
                        ref_mh.add_hash(hash);
                    }
                }
                // Calculating similarity and AAI
                let similarity = query_mh.similarity(&ref_mh, false, false).unwrap_or(0.0);
                let aai = 1.0 + ((2.0 * similarity) / (1.0 + similarity)).ln() / kmer_size as f64;
                format!("{}\t{}\t{}\t{}", query_file, ref_file, similarity, aai)
            })
        })
        .collect();

    std::fs::write(output_file_path, results.join("\n"))?;
    Ok(())
}
