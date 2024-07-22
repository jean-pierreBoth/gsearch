use clap::{Command, Arg, value_parser};
use hmmer_rs::{Hmm, HmmerPipeline, EaselSequence, Alphabet};
use needletail::parse_fastx_file;
use std::path::Path;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::ffi::CStr;
use std::error::Error;
use env_logger::Builder;
use std::convert::TryInto; // Make sure to include TryInto

pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    1
}

fn main() -> Result<(), Box<dyn Error>> {
    init_log();
    let start_t = chrono::Local::now();
    log::info!("\n hmmsearch begins at time:{:#?} \n ", start_t);

    let matches = Command::new("hmmsearch_rs")
        .arg_required_else_help(true)
        .about("Search protein sequences against HMM profiles")
        .version("0.1.0")
        .arg(Arg::new("fasta")
            .short('f')
            .long("faa")
            .help("Path to the faa file containing the protein sequences")
            .required(true)
            .value_parser(value_parser!(String)))
        .arg(Arg::new("hmm")
            .short('m')
            .long("hmm")
            .help("Path to the HMM file")
            .required(true)
            .value_parser(value_parser!(String)))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .help("Output file to save the search results")
            .required(false)
            .value_parser(value_parser!(String)))
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap().clone();
    let hmm_path = matches.get_one::<String>("hmm").unwrap().clone();
    let output_path = matches.get_one::<String>("output")
        .map_or_else(|| "output.txt".to_string(), |s| s.clone());

    let hmms = Hmm::read_hmms_from_path(Path::new(&hmm_path))?;

    let mut reader = parse_fastx_file(Path::new(&fasta_path))?;
    let mut results: Vec<(f32, String, usize, String, f32, f64, String)> = Vec::new();

    while let Some(Ok(record)) = reader.next() {
        let seq_id = String::from_utf8(record.id().to_vec())?;
        let sequence = record.seq().to_vec();

        for hmm in &hmms {
            let mut hmmsearch = HmmerPipeline::new(hmm);
            let mut query_seq = EaselSequence::new(Alphabet::Protein);
            query_seq.replace_sequence(&sequence)?;
            hmmsearch.query(&query_seq);
            let hmmsearch_result = hmmsearch.get_results();
            drop(hmmsearch);
            unsafe {
                if !hmmsearch_result.c_th.is_null() {
                    let th = &*hmmsearch_result.c_th;
                    if th.nreported > 0 && !th.hit.is_null() {
                        let first_hit = &**th.hit.offset(0);
                        if !(*first_hit).dcl.is_null() {
                            let first_hit_name = CStr::from_ptr((*first_hit).name).to_string_lossy();
                            let first_hit_score = (*first_hit).score;
                            let first_domain = &*(*first_hit).dcl.offset(0);
                            let first_domain_score = first_domain.bitscore;
                            let first_domain_evalue = first_domain.lnP.exp() * (*hmmsearch_result.c_pli).Z;

                            results.push((
                                first_hit_score,
                                hmm.name().to_string(),
                                th.nreported.try_into().unwrap(), // Convert here
                                first_hit_name.to_string(),
                                first_domain_score,
                                first_domain_evalue,
                                seq_id.clone()
                            ));
                        }
                    }
                }
            }
        }
    }

    let mut output_file = BufWriter::new(File::create(&output_path)?);
    writeln!(output_file, "HMM\tScore\tNreported\tHit Name\tDomain Score\tE-value\tSequence")?;
    for (score, hmm_name, nreported, first_hit_name, domain_score, domain_evalue, seq_id) in &results {
        writeln!(
            output_file,
            "{}\t{:.4}\t{}\t{}\t{:.4}\t{:.2e}\t{}",
            hmm_name, score, nreported, first_hit_name, domain_score, domain_evalue, seq_id
        )?;
    }

    let end_t = chrono::Local::now();
    log::info!("\n hmmsearch ends at time:{:#?} \n ", end_t);

    Ok(())
}
