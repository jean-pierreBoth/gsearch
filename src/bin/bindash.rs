use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::collections::HashMap;
use std::fmt::Debug;
use std::path::Path;

use needletail::{parse_fastx_file, Sequence};
use num;

use kmerutils::sketcharg::{SeqSketcherParams, SketchAlgo, DataType};
use kmerutils::base::{
    alphabet::Alphabet2b,
    sequence::Sequence as SequenceStruct,
    kmergenerator::*,
    Kmer32bit,
    Kmer16b32bit,
    Kmer64bit,
    CompressedKmerT,
    KmerBuilder,
};
use kmerutils::sketching::setsketchert::{
    SeqSketcherT, // trait
    OptDensHashSketch,
    RevOptDensHashSketch,
};

use anndists::dist::{Distance, DistHamming};

/// Converts ASCII-encoded bases (from Needletail) into our `SequenceStruct`.
fn ascii_to_seq(bases: &[u8]) -> Result<SequenceStruct, ()> {
    let alphabet = Alphabet2b::new();
    let mut seq = SequenceStruct::with_capacity(2, bases.len());
    seq.encode_and_add(bases, &alphabet);
    Ok(seq)
}

/// Reads a list of file paths (one per line) from a text file.
fn read_genome_list(filepath: &str) -> Vec<String> {
    let file = File::open(filepath).expect("Cannot open genome list file");
    let reader = BufReader::new(file);
    reader
        .lines()
        .map(|line| line.expect("Error reading genome list"))
        .collect()
}

/// Generic function that sketches a list of FASTA/Q file paths using a provided
/// `SeqSketcherT` and a k-mer hash function. Returns a `HashMap` from file path
/// to the single sketch vector.
///
/// ### **Key**: The `sketcher` reference now has `+ Sync` so it can be shared
/// across threads. We also pass `kmer_hash_fn` by value (requires `Copy + Send + Sync`).
fn sketch_files<Kmer, F>(
    file_paths: &[String],
    sketcher: &(impl SeqSketcherT<Kmer, Sig = f32> + Sync), 
    kmer_hash_fn: F,
) -> HashMap<String, Vec<f32>>
where
    // Required for using `SeqSketcherT` with Kmer
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,

    // The function or closure must be Copy + Send + Sync to work in parallel
    F: Fn(&Kmer) -> <Kmer as CompressedKmerT>::Val + Send + Sync + Copy,
{
    file_paths
        .par_iter()
        .map(|path| {
            // Read all sequences from this FASTA/Q
            let mut sequences = Vec::new();
            let mut reader = parse_fastx_file(path).expect("Invalid FASTA/Q file");
            while let Some(record) = reader.next() {
                let seq_record = record.expect("Error reading sequence record");
                let seq_seq = seq_record.normalize(false).into_owned();
                let seq = ascii_to_seq(&seq_seq).unwrap();
                sequences.push(seq);
            }

            // Prepare references
            let sequences_ref: Vec<&SequenceStruct> = sequences.iter().collect();

            // We call the sketcher method
            let signature = sketcher.sketch_compressedkmer_seqs(&sequences_ref, kmer_hash_fn);

            // Return the single signature for the file
            (path.clone(), signature[0].clone())
        })
        .collect::<HashMap<String, Vec<f32>>>()
}

/// Computes the distance between two sketches (as `Vec<f32>`) based on
/// Hamming in the densified sketches, then uses your transformation:
///
///    `distance = -ln( (2*j) / (1 + j) ) / kmer_size`.
fn compute_distance(query_sig: &[f32], reference_sig: &[f32], kmer_size: usize) -> f64 {
    let dist_hamming = DistHamming;
    let hamming_distance = dist_hamming.eval(query_sig, reference_sig);
    let hamming_distance = if hamming_distance == 0.0 {
        std::f32::EPSILON // Use a small value close to zero
    } else {
        hamming_distance
    };
    // Jaccard from Hamming
    let j = 1.0 - hamming_distance;
    let numerator = 2.0 * j;
    let denominator = 1.0 + j;
    let fraction = (numerator as f64) / (denominator as f64);

    -fraction.ln() / (kmer_size as f64)
}

fn write_results(
    output: Option<String>,
    query_genomes: &[String],
    reference_genomes: &[String],
    query_sketches: &HashMap<String, Vec<f32>>,
    reference_sketches: &HashMap<String, Vec<f32>>,
    kmer_size: usize,
) {
    let mut output_writer: Box<dyn Write> = match output {
        Some(filename) => {
            Box::new(BufWriter::new(File::create(&filename).expect("Cannot create output file")))
        }
        None => Box::new(BufWriter::new(io::stdout())),
    };

    writeln!(output_writer, "Query\tReference\tDistance").expect("Error writing header");

    // 1) Build list of all query/reference pairs
    //    (Potentially large, but straightforward to parallelize).
    let all_pairs: Vec<(String, String)> = query_genomes
        .iter()
        .flat_map(|q| {
            reference_genomes.iter().map(move |r| (q.clone(), r.clone()))
        })
        .collect();

    // 2) Parallel compute over those pairs
    let results: Vec<(String, String, f64)> = all_pairs
        .into_par_iter()
        .map(|(query_path, reference_path)| {
            let query_signature = &query_sketches[&query_path];
            let reference_signature = &reference_sketches[&reference_path];

            let distance = compute_distance(query_signature, reference_signature, kmer_size);

            // If the file basename is the same, consider distance=0.0
            let query_basename = Path::new(&query_path)
                .file_name()
                .and_then(|os_str| os_str.to_str())
                .unwrap_or(&query_path);

            let reference_basename = Path::new(&reference_path)
                .file_name()
                .and_then(|os_str| os_str.to_str())
                .unwrap_or(&reference_path);

            let distance = if query_basename == reference_basename {
                0.0
            } else {
                distance
            };

            (query_path, reference_path, distance)
        })
        .collect();

    // 3) Write out the results (in whatever order they were computed)
    for (q_path, r_path, dist) in results {
        writeln!(output_writer, "{}\t{}\t{:.6}", q_path, r_path, dist)
            .expect("Error writing output");
    }
}

/// Runs the pipeline for a given Kmer type and densification type (`dens`).
///
/// - `dens = 0` => `OptDensHashSketch`
/// - `dens = 1` => `RevOptDensHashSketch`
///
/// The `kmer_hash_fn` is taken by value (and must be `Copy + Send + Sync`).
fn sketching_kmerType<Kmer, F>(
    query_genomes: &[String],
    reference_genomes: &[String],
    sketch_args: &SeqSketcherParams,
    kmer_hash_fn: F,
    dens: usize,
    output: Option<String>,
    kmer_size: usize,
) where
    // The Kmer must meet the constraints for the chosen sketchers
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,

    // The hash function must be usable in parallel
    F: Fn(&Kmer) -> <Kmer as CompressedKmerT>::Val + Send + Sync + Copy,
{
    match dens {
        0 => {
            // dens=0 => "optimal" densification
            let sketcher = OptDensHashSketch::<Kmer, f32>::new(sketch_args);
            println!("Sketching query genomes with Optimal Densification MinHash...");
            let query_sketches = sketch_files(query_genomes, &sketcher, kmer_hash_fn);

            println!("Sketching reference genomes Optimal Densification MinHash...");
            let reference_sketches = sketch_files(reference_genomes, &sketcher, kmer_hash_fn);

            println!("Performing pairwise comparisons...");
            write_results(
                output,
                query_genomes,
                reference_genomes,
                &query_sketches,
                &reference_sketches,
                kmer_size,
            );
        }
        1 => {
            // dens=1 => "reverse optimal" densification
            let sketcher = RevOptDensHashSketch::<Kmer, f32>::new(sketch_args);
            println!("Sketching query genomes with Reverse Optimal Densification MinHash...");
            let query_sketches = sketch_files(query_genomes, &sketcher, kmer_hash_fn);

            println!("Sketching reference genomes with Revse Optimal Densification MinHash...");
            let reference_sketches = sketch_files(reference_genomes, &sketcher, kmer_hash_fn);

            println!("Performing pairwise comparisons...");
            write_results(
                output,
                query_genomes,
                reference_genomes,
                &query_sketches,
                &reference_sketches,
                kmer_size,
            );
        }
        _ => {
            panic!("Only densification = 0 or 1 are supported!");
        }
    }
}

/// Main entry point of the program.
fn main() {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();

    let matches = Command::new("BinDash")
        .version("0.2.9")
        .about("Binwise Densified MinHash for Genome/Metagenome/Pangenome Comparisons")
        .arg(
            Arg::new("query_list")
                .short('q')
                .long("query_list")
                .value_name("QUERY_LIST_FILE")
                .help("Query genome list file (one FASTA/FNA file path per line, .gz supported)")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("reference_list")
                .short('r')
                .long("reference_list")
                .value_name("REFERENCE_LIST_FILE")
                .help("Reference genome list file (one FASTA/FNA file path per line, .gz supported)")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("kmer_size")
                .short('k')
                .long("kmer_size")
                .value_name("KMER_SIZE")
                .help("K-mer size")
                .default_value("16")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("sketch_size")
                .short('s')
                .long("sketch_size")
                .value_name("SKETCH_SIZE")
                .help("MinHash sketch size")
                .default_value("2048")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("dens_opt")
                .short('d')
                .long("densification")
                .value_name("DENS_OPT")
                .help("Densification strategy, 0 = optimal densification, 1 = reverse optimal/faster densification")
                .default_value("0")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads to use in parallel")
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_FILE")
                .help("Output file (defaults to stdout)")
                .required(false)
                .action(ArgAction::Set),
        )
        .get_matches();

    let query_list = matches
        .get_one::<String>("query_list")
        .expect("Query list is required")
        .to_string();
    let reference_list = matches
        .get_one::<String>("reference_list")
        .expect("Reference list is required")
        .to_string();
    let kmer_size = *matches.get_one::<usize>("kmer_size").unwrap();
    let sketch_size = *matches.get_one::<usize>("sketch_size").unwrap();
    let dens = *matches.get_one::<usize>("dens_opt").unwrap();
    let threads = *matches.get_one::<usize>("threads").unwrap();
    let output = matches.get_one::<String>("output").cloned();

    // Build global threadpool
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    // Read genome lists
    let query_genomes = read_genome_list(&query_list);
    let reference_genomes = read_genome_list(&reference_list);

    // Build the base sketch arguments
    // We'll pick the densification type via `dens` below.
    let sketch_args = SeqSketcherParams::new(
        kmer_size,
        sketch_size,
        SketchAlgo::OPTDENS, // placeholder
        DataType::DNA,
    );

    // Depending on kmer_size, pick a Kmer type and define the hash function (by value).
    if kmer_size <= 14 {
        // Kmer32bit
        let nb_alphabet_bits = 2;
        let kmer_hash_fn_32bit = move |kmer: &Kmer32bit| -> <Kmer32bit as CompressedKmerT>::Val {
            let mask: <Kmer32bit as CompressedKmerT>::Val =
                num::NumCast::from::<u64>((1u64 << (nb_alphabet_bits * kmer.get_nb_base())) - 1)
                    .unwrap();
            kmer.get_compressed_value() & mask
        };

        sketching_kmerType::<Kmer32bit, _>(
            &query_genomes,
            &reference_genomes,
            &sketch_args,
            kmer_hash_fn_32bit,
            dens,
            output,
            kmer_size,
        );

    } else if kmer_size == 16 {
        // Kmer16b32bit
        let nb_alphabet_bits = 2;
        let kmer_hash_fn_16b32bit = move |kmer: &Kmer16b32bit| -> <Kmer16b32bit as CompressedKmerT>::Val {
            // canonical
            let canonical = kmer.reverse_complement().min(*kmer);
            let mask: <Kmer16b32bit as CompressedKmerT>::Val =
                num::NumCast::from::<u64>((1u64 << (nb_alphabet_bits * kmer.get_nb_base())) - 1)
                    .unwrap();
            canonical.get_compressed_value() & mask
        };

        sketching_kmerType::<Kmer16b32bit, _>(
            &query_genomes,
            &reference_genomes,
            &sketch_args,
            kmer_hash_fn_16b32bit,
            dens,
            output,
            kmer_size,
        );

    } else if kmer_size <= 32 {
        // Kmer64bit
        let nb_alphabet_bits = 2;
        let kmer_hash_fn_64bit = move |kmer: &Kmer64bit| -> <Kmer64bit as CompressedKmerT>::Val {
            let canonical = kmer.reverse_complement().min(*kmer);
            let mask: <Kmer64bit as CompressedKmerT>::Val =
                num::NumCast::from::<u64>((1u64 << (nb_alphabet_bits * kmer.get_nb_base())) - 1)
                    .unwrap();
            canonical.get_compressed_value() & mask
        };

        sketching_kmerType::<Kmer64bit, _>(
            &query_genomes,
            &reference_genomes,
            &sketch_args,
            kmer_hash_fn_64bit,
            dens,
            output,
            kmer_size,
        );

    } else {
        panic!("kmer_size must not be 15 and cannot exceed 32!");
    }
}
