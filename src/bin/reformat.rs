use clap::{Arg, Command};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

fn main() -> io::Result<()> {
    let matches = Command::new("reformat")
        .version("0.2.9")
        .author("Your Name")
        .about("Processes input files for ANI calculation")
        .arg(Arg::new("kmer")
             .help("The kmer value used for ANI calculation")
             .required(true))
        .arg(Arg::new("model")
             .help("The model to be used for ANI calculation (1 or 2, corresponding to Poisson model and Binomial model)")
             .required(true))
        .arg(Arg::new("input_file")
             .help("File containing the data to be transformed into tabular format")
             .required(true))
        .arg(Arg::new("output_file")
             .help("File where the tabular output will be saved")
             .required(true))
        .get_matches();

    let kmer: i32 = matches.get_one::<String>("kmer")
                           .unwrap()
                           .parse::<i32>()
                           .expect("Kmer must be an integer");
    let model: i32 = matches.get_one::<String>("model")
                           .unwrap()
                           .parse::<i32>()
                           .expect("Model must be an integer");
    let input_file = matches.get_one::<String>("input_file").unwrap().as_str();
    let output_file = matches.get_one::<String>("output_file").unwrap().as_str();
    
    // Read and process the input file
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);

    let mut results: Vec<String> = reader
        .lines()
        .filter_map(Result::ok)
        .par_bridge()
        .filter(|line| line.starts_with("query_id:"))
        .map(|line| process_line(&line, kmer, model))
        .collect();

    // Parallel sort the results, excluding the header
    results.par_sort_by(|a, b| {
        let a_cols: Vec<&str> = a.split('\t').collect();
        let b_cols: Vec<&str> = b.split('\t').collect();
        a_cols[0].cmp(b_cols[0])
            .then_with(|| {
                a_cols[1].parse::<f64>().unwrap_or_default().partial_cmp(&b_cols[1].parse::<f64>().unwrap_or_default()).unwrap()
            })
    });

    // Write output
    let mut output = File::create(output_file)?;
    writeln!(output, "Query_Name\tDistance\tNeighbor_Fasta_name\tNeighbor_Seq_Len\tANI")?;
    for line in results {
        writeln!(output, "{}", line)?;
    }

    Ok(())
}

fn process_line(line: &str, kmer: i32, model: i32) -> String {
    let parts: Vec<&str> = line.split('\t').collect();
    let query_id = Path::new(parts[1]).file_name().unwrap().to_str().unwrap().to_string();
    let distance_str = parts[3];
    let distance: f64 = distance_str.parse().unwrap_or_else(|_| panic!("Failed to parse distance '{}' in line: {}", distance_str, line));
    let answer_fasta_path = Path::new(parts[5]).file_name().unwrap().to_str().unwrap().to_string();
    let answer_seq_len = parts[7];
    let ani = calculate_ani(distance, kmer, model);
    format!("{}\t{}\t{}\t{}\t{}", query_id, distance, answer_fasta_path, answer_seq_len, ani)
}

fn calculate_ani(distance: f64, kmer: i32, model: i32) -> String {
    match model {
        1 => format!("{}", (1.0 + ((1.0 - distance) * 2.0 / (1.0 - distance + 1.0)).ln() / kmer as f64) * 100.0),
        2 => format!("{}", ((1.0 - distance) * 2.0 / (1.0 - distance + 1.0)).powf(1.0 / kmer as f64) * 100.0),
        _ => "Invalid Model".to_string(),
    }
}
