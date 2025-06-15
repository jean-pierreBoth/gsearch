use clap::{Arg, Command};
use bigsig::{bigsi, build};
use std::alloc::System;
use std::io::Result;
use std::time::SystemTime;
use env_logger::Builder;
use rayon::ThreadPoolBuilder;

pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    1
}

#[global_allocator]
static GLOBAL: System = System;

fn main() -> Result<()> {
    let _ = init_log();

    let matches = Command::new("bigsig")
        .version("0.2.9")
        .about("Large-scale Sequence Search with BItsliced Genomic Signature Index (BIGSIG)")
        .arg_required_else_help(true)
        .subcommand(Command::new("construct")
            .about("Construct a BIGSIG")
            .version("0.1.0")
            .arg_required_else_help(true)
            .arg(Arg::new("bigsi")
                .short('b')
                .long("bigsi")
                .value_name("FILE")
                .help("BIGSI file to output")
                .required(true))
            .arg(Arg::new("ref_file")
                .short('r')
                .long("refs")
                .value_name("FILE")
                .help("Sets the input reference file to use")
                .required(true))
            .arg(Arg::new("k-mer_size")
                .short('k')
                .long("kmer")
                .value_name("INT")
                .help("Sets k-mer size")
                .required(true))
            .arg(Arg::new("num_hashes")
                .short('n')
                .long("num_hashes")
                .value_name("INT")
                .help("Sets number of hashes for bloom filter")
                .required(true))
            .arg(Arg::new("length_bloom")
                .short('s')
                .long("bloom")
                .value_name("INT")
                .help("Sets the length of the bloom filter")
                .required(true))
            .arg(Arg::new("minimizer")
                .short('m')
                .long("minimizer")
                .action(clap::ArgAction::SetTrue)
                .help("Build index with minimizers"))
            .arg(Arg::new("value")
                .short('v')
                .long("value")
                .value_name("INT")
                .help("Sets length minimizer (default 15)")
                .default_value("15"))
            .arg(Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("INT")
                .help("Number of threads to use")
                .default_missing_value("1"))
            .arg(Arg::new("quality")
                .short('Q')
                .long("quality")
                .value_name("INT")
                .help("Minimum phred score to keep basepairs within read")
                .default_value("15"))
            .arg(Arg::new("filter")
                .short('f')
                .long("filter")
                .value_name("INT")
                .help("Minimum coverage kmer threshold")
                .default_missing_value("-1")))
        .subcommand(Command::new("query")
            .about("Query a BIGSIG on one or more fasta/fastq.gz files")
            .arg_required_else_help(true)
            .version("0.1.0")
            .arg(Arg::new("bigsi")
                .short('b')
                .long("bigsi")
                .value_name("FILE")
                .help("Sets the name of the index file for search")
                .required(true))
            .arg(Arg::new("query")
                .short('q')
                .long("query")
                .help("Query file(-s) fastq.gz")
                .required(true)
                .num_args(1..))
            .arg(Arg::new("reverse")
                .short('r')
                .long("reverse")
                .value_name("FILE")
                .help("Reverse file(-s) fastq.gz")
                .required(false)
                .num_args(1..)
                .default_value("none"))
            .arg(Arg::new("filter")
                .short('f')
                .long("filter")
                .value_name("INT")
                .help("Set minimum k-mer frequency")
                .required(false))
            .arg(Arg::new("shared_kmers")
                .short('p')
                .long("p_shared")
                .value_name("FLOAT")
                .help("Set minimum proportion of shared k-mers with a reference")
                .required(false))
            .arg(Arg::new("gene_search")
                .short('g')
                .long("gene_search")
                .action(clap::ArgAction::SetTrue)
                .help("If set, the proportion of kmers from the query matching the entries in the index will be reported"))
            .arg(Arg::new("perfect_search")
                .short('s')
                .long("perfect_search")
                .action(clap::ArgAction::SetTrue)
                .help("If set, the fast 'perfect match' algorithm will be used"))
            .arg(Arg::new("multi_fasta")
                .short('m')
                .long("multi_fasta")
                .action(clap::ArgAction::SetTrue)
                .help("If set, each accession in a multifasta will be treated as a separate query, currently only with the -s option"))
            .arg(Arg::new("quality")
                .short('Q')
                .long("quality")
                .value_name("INT")
                .help("Minimum phred score to keep basepairs within read (default 15)")
                .default_value("15")))
        .subcommand(Command::new("identify")
            .about("Identify reads based on probability")
            .version("0.1.0")
            .arg(Arg::new("bigsi")
                .short('b')
                .long("bigsi")
                .help("Index to be used for search")
                .required(true))
            .arg(Arg::new("query")
                .short('q')
                .long("query")
                .help("Query file(-s) fastq.gz")
                .required(true)
                .num_args(1..)) 
            .arg(Arg::new("batch")
                .short('c')
                .long("batch")
                .help("Sets size of batch of reads to be processed in parallel (default 50,000)")
                .required(false)
                .default_value("50000"))
            .arg(Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("INT")
                .help("Number of threads to use, if not set the maximum available number threads will be used")
                .required(false)
                .default_value("0"))
            .arg(Arg::new("prefix")
                .short('n')
                .long("prefix")
                .help("Prefix for output file(-s)")
                .required(true))
            .arg(Arg::new("down_sample")
                .short('d')
                .long("down_sample")
                .help("Down-sample k-mers used for read classification, default 1; increases speed at cost of decreased sensitivity")
                .default_value("1"))
            .arg(Arg::new("high_mem_load")
                .short('H')
                .long("high_mem_load")
                .action(clap::ArgAction::SetTrue)
                .help("When this flag is set, a faster, but less memory efficient method to load the index is used"))
            .arg(Arg::new("fp_correct")
                .short('p')
                .long("fp_correct")
                .required(false)
                .help("Parameter to correct for false positives, default 3 (= 0.001), may be increased for larger searches")
                .default_value("3.0"))
            .arg(Arg::new("quality")
                .short('Q')
                .long("quality")
                .help("kmers with nucleotides below this minimum phred score will be excluded from the analyses (default 15)")
                .default_value("15"))
            .arg(Arg::new("bitvector_sample")
                .short('B')
                .long("bitvector_sample")
                .help("Collects matches for subset of kmers indicated (default=3), using this subset to more rapidly find hits for the remainder of the kmers")
                .default_value("3")))
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("construct") {
        let ref_file = matches.get_one::<String>("ref_file").expect("required");
        let bigsi_file = matches.get_one::<String>("bigsi").expect("required");
        let kmer_size = matches.get_one::<String>("k-mer_size").expect("required").parse::<usize>().expect("Integer required");
        let num_hashes = matches.get_one::<String>("num_hashes").expect("required").parse::<usize>().expect("Integer required");
        let length_bloom = matches.get_one::<String>("length_bloom").expect("required").parse::<usize>().expect("Integer required");
        let threads = matches.get_one::<String>("threads").unwrap_or(&"0".to_string()).parse::<usize>().expect("Integer required");
        let quality = matches.get_one::<String>("quality").unwrap_or(&"15".to_string()).parse::<u8>().expect("Integer required");
        let minimizer = matches.contains_id("minimizer");
        let minimizer_value = matches.get_one::<String>("value").unwrap_or(&"21".to_string()).parse::<usize>().expect("Integer required");
        let filter = matches.get_one::<String>("filter").unwrap_or(&"-1".to_string()).parse::<isize>().expect("Integer required");

        let map = build::tab_to_map(ref_file.to_string());
        
        if minimizer {
            println!("Building with minimizers, minimizer size: {}", minimizer_value);
            let (bigsi_map, colors_accession, n_ref_kmers) = if threads == 1 {
                build::build_single_mini(&map, length_bloom, num_hashes, kmer_size, minimizer_value, quality, filter)
            } else {
                build::build_multi_mini(&map, length_bloom, num_hashes, kmer_size, minimizer_value, threads, quality, filter)
            };
            println!("Saving BIGSI to file: {}", bigsi_file);
            bigsi::save_bigsi_mini(bigsi_file, &bigsi::BigsyMapMiniNew {
                map: bigsi_map,
                colors: colors_accession,
                n_ref_kmers: n_ref_kmers,
                bloom_size: length_bloom,
                num_hash: num_hashes,
                k_size: kmer_size,
                m_size: minimizer_value,
            });
        } else {
            let (bigsi_map, colors_accession, n_ref_kmers) = if threads == 1 {
                build::build_single(&map, length_bloom, num_hashes, kmer_size, quality, filter)
            } else {
                build::build_multi(&map, length_bloom, num_hashes, kmer_size, threads, quality, filter)
            };
            println!("Saving BIGSI to file: {}", bigsi_file);
            bigsi::save_bigsi(bigsi_file, &bigsi::BigsyMapNew {
                map: bigsi_map,
                colors: colors_accession,
                n_ref_kmers: n_ref_kmers,
                bloom_size: length_bloom,
                num_hash: num_hashes,
                k_size: kmer_size,
            });
        }
    }
    if let Some(matches) = matches.subcommand_matches("query") {
        let files1: Vec<&str> = matches.get_many::<String>("query").unwrap().map(|s| s.as_str()).collect();
        let files2: Vec<&str> = if matches.get_one::<String>("reverse").unwrap() == "none" {
            vec![]
        } else {
            matches.get_many::<String>("reverse").unwrap().map(|s| s.as_str()).collect()
        };
        let filter = matches.get_one::<String>("filter").unwrap_or(&"-1".to_string()).parse::<isize>().unwrap();
        let cov = matches.get_one::<String>("shared_kmers").unwrap_or(&"0.35".to_string()).parse::<f64>().unwrap();
        let gene_search = matches.contains_id("gene_search");
        let perfect_search = matches.contains_id("perfect_search");
        let multi_fasta = matches.contains_id("multi_fasta");
        let quality = matches.get_one::<String>("quality").unwrap_or(&"15".to_string()).parse::<u8>().unwrap();
        if matches.get_one::<String>("bigsi").unwrap().ends_with(".mxi") {
            eprintln!(
                "Error: An index with minimizers (.mxi) is used, but not available for this function"
            );
        } else {
            let bigsi_time = SystemTime::now();
            eprintln!("Loading index");
            let index = bigsi::read_bigsi(matches.get_one::<String>("bigsi").unwrap());
            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if perfect_search {
                if multi_fasta {
                    //make 'perfect' batch....
                    bigsig::perfect_search::batch_search_mf(
                        files1,
                        &index.map,
                        &index.colors,
                        &index.n_ref_kmers,
                        index.bloom_size,
                        index.num_hash,
                        index.k_size,
                        cov,
                    )
                } else {
                    bigsig::perfect_search::batch_search(
                        files1,
                        &index.map,
                        &index.colors,
                        &index.n_ref_kmers,
                        index.bloom_size,
                        index.num_hash,
                        index.k_size,
                        cov,
                    )
                }
            } else {
                bigsig::batch_search_pe::batch_search(
                    files1,
                    files2,
                    &index.map,
                    &index.colors,
                    &index.n_ref_kmers,
                    index.bloom_size,
                    index.num_hash,
                    index.k_size,
                    filter,
                    cov,
                    gene_search,
                    quality,
                )
            }
        }
    }
    if let Some(matches) = matches.subcommand_matches("identify") {
        let bigsi_time = SystemTime::now();
        let fq: Vec<&str> = matches.get_many::<String>("query").unwrap().map(|s| s.as_str()).collect();
        let threads: usize = matches.get_one::<String>("threads").unwrap_or(&"0".to_string()).parse().expect("Invalid threads number");
        let down_sample: usize = matches.get_one::<String>("down_sample")
            .map(|s| s.parse::<usize>().expect("Invalid down_sample value"))
            .unwrap_or(1);
        let correct: f64 = matches.get_one::<String>("fp_correct").unwrap_or(&"3.0".to_string()).parse().expect("Invalid fp_correct value");
        let fp_correct = 10f64.powf(-correct);
        let index = matches.get_one::<String>("bigsi").unwrap();
        let prefix = matches.get_one::<String>("prefix").unwrap();
        let quality: u8 = matches.get_one::<String>("quality").unwrap_or(&"15".to_string()).parse().expect("Invalid quality value");
        let batch: usize = matches.get_one::<String>("batch").unwrap_or(&"50000".to_string()).parse().expect("Invalid batch value");
        let high_mem_load = matches.contains_id("high_mem_load");
        let bitvector_sample: usize = matches.get_one::<String>("bitvector_sample").unwrap_or(&"3".to_string()).parse().expect("Invalid bitvector_sample value");
        
        ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .expect("Can't initialize ThreadPoolBuilder");

        if index.ends_with(".mxi") {
            //let metadata = fs::metadata(&index).expect("Can't read metadata index!");
            let bigsi = if high_mem_load {
                bigsi::read_bigsi_mini_highmem(index)
            } else {
                bigsi::read_bigsi_mini(index)
            };
            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if fq[0].ends_with(".gz") {
                if fq.len() > 1 {
                    bigsig::read_id_mt_pe::per_read_stream_pe(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        bigsi.m_size,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                        bitvector_sample,
                    )
                } else {
                    bigsig::read_id_mt_pe::per_read_stream_se(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        bigsi.m_size,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                        bitvector_sample,
                    )
                };
            } else {
                bigsig::read_id_mt_pe::stream_fasta(
                    fq,
                    &bigsi.map,
                    &bigsi.colors,
                    &bigsi.n_ref_kmers,
                    bigsi.bloom_size,
                    bigsi.num_hash,
                    bigsi.k_size,
                    bigsi.m_size,
                    threads,
                    down_sample,
                    fp_correct,
                    batch,
                    prefix,
                    bitvector_sample,
                );
            }
        } else {
            let bigsi = if high_mem_load {
                bigsi::read_bigsi_highmem(index)
            } else {
                bigsi::read_bigsi(index)
            };
            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if fq[0].ends_with(".gz") {
                if fq.len() > 1 {
                    bigsig::read_id_mt_pe::per_read_stream_pe(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        0,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                        bitvector_sample,
                    )
                } else {
                    bigsig::read_id_mt_pe::per_read_stream_se(
                        fq,
                        &bigsi.map,
                        &bigsi.colors,
                        &bigsi.n_ref_kmers,
                        bigsi.bloom_size,
                        bigsi.num_hash,
                        bigsi.k_size,
                        0,
                        threads,
                        down_sample,
                        fp_correct,
                        batch,
                        prefix,
                        quality,
                        bitvector_sample,
                    )
                };
            } else {
                bigsig::read_id_mt_pe::stream_fasta(
                    fq,
                    &bigsi.map,
                    &bigsi.colors,
                    &bigsi.n_ref_kmers,
                    bigsi.bloom_size,
                    bigsi.num_hash,
                    bigsi.k_size,
                    0,
                    threads,
                    down_sample,
                    fp_correct,
                    batch,
                    prefix,
                    bitvector_sample,
                );
            }
        }
        bigsig::reports::read_counts_five_fields(prefix.to_owned() + "_reads.txt", prefix);
    }
    Ok(())
}
