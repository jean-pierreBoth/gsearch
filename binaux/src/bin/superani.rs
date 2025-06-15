use clap::{Arg, Command};
use env_logger::Builder;
use rayon::prelude::*;
use skani::chain;
use skani::file_io;
use skani::params::*;
use skani::regression;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::sync::{Arc, Mutex};

fn default_params(mode: Mode) -> (Arc<CommandParams>, Arc<SketchParams>) {
    let cmd_params = Arc::new(CommandParams {
        screen: true,
        screen_val: 75.00,
        mode,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        max_results: 10000000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: 0.10,
        keep_refs: false,
        est_ci: false,
        learned_ani: true,
        detailed_out: false,
        diagonal: false,
        rescue_small: false,
        distance: false,
    });

    let m = 1000;
    let c = 30;
    let k = 16;
    let sketch_params = Arc::new(SketchParams::new(m, c, k, false, false));
    (cmd_params, sketch_params)
}

pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    1
}

fn main() -> io::Result<()> {
    let _ = init_log();
    let start_t = chrono::Local::now();
    log::info!("\n superani begins at time:{:#?} \n ", start_t);

    let matches = Command::new("SuperANI")
    .about("Computing average nucleotide identity between reference and query genomes via sparse kmer chaining or Open Syncmer with Densified MinHash")
    .version("0.2.9")
    .arg(
        Arg::new("query_list")
            .short('q')
            .long("ql")
            .value_name("FILE")
            .help("A file containing a list of query genome paths (.gz supported)")
            .required(true)
            .value_parser(clap::value_parser!(String))
    )
    .arg(
        Arg::new("ref_list")
            .short('r')
            .long("rl")
            .value_name("FILE")
            .help("A file containing a list of reference genome paths (.gz supported)")
            .required(true)
            .value_parser(clap::value_parser!(String))
    )
    .arg(
        Arg::new("output")
            .short('o')
            .long("output")
            .value_name("FILE")
            .help("Output file to write results")
            .required(true)
            .value_parser(clap::value_parser!(String))
    )
    .get_matches();

    let query_file_path = matches.get_one::<String>("query_list").unwrap();
    let ref_file_path = matches.get_one::<String>("ref_list").unwrap();
    let output_file_path = matches.get_one::<String>("output").unwrap();

    let queries = Arc::new(
        BufReader::new(File::open(query_file_path)?)
            .lines()
            .collect::<Result<Vec<_>, io::Error>>()?,
    );
    let refs = Arc::new(
        BufReader::new(File::open(ref_file_path)?)
            .lines()
            .collect::<Result<Vec<_>, io::Error>>()?,
    );

    let output_file = Arc::new(Mutex::new(File::create(output_file_path)?));

    let (command_params, sketch_params) = default_params(Mode::Dist);
    let model_opt = regression::get_model(sketch_params.c, true);

    refs.par_iter().for_each(|ref_file| {
        let ref_sketches =
            file_io::fastx_to_sketches(&vec![ref_file.clone()], &sketch_params, true);
        let queries_clone = queries.clone();
        let output_file_clone = output_file.clone();
        let command_params_clone = command_params.clone();
        // let model_opt_clone = model_opt.clone();
        let sketch_params_clone = sketch_params.clone(); // Clone Arc for use inside the nested closure
        queries_clone.par_iter().for_each(|query_file| {
            let query_sketches =
                file_io::fastx_to_sketches(&vec![query_file.clone()], &sketch_params_clone, true);
            for ref_sketch in ref_sketches.iter() {
                for query_sketch in query_sketches.iter() {
                    let map_params = chain::map_params_from_sketch(
                        ref_sketch,
                        false,
                        &command_params_clone,
                        &model_opt,
                    );
                    let mut ani_result = chain::chain_seeds(ref_sketch, query_sketch, map_params);
                    if let Some(model) = model_opt.as_ref() {
                        regression::predict_from_ani_res(&mut ani_result, model);
                    }
                    let output = format!(
                        "{}\t{}\t{}\t{}\t{}\n",
                        query_file,
                        ref_file,
                        ani_result.ani,
                        ani_result.align_fraction_query,
                        ani_result.align_fraction_ref
                    );
                    let mut file = output_file_clone.lock().unwrap();
                    file.write_all(output.as_bytes()).unwrap();
                }
            }
        });
    });
    Ok(())
}
