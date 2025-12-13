use clap::{Arg, ArgAction, Command};
use gsearch::utils::reloadhnsw;
use gsearch::utils::SeqDict;
use log::info;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use annembed::fromhnsw::kgraph::KGraph;
use annembed::fromhnsw::kgraph_from_hnsw_all;
use hnsw_rs::prelude::DistHamming;
use kmerutils::base::Kmer32bit;
use kmerutils::sketching::setsketchert::*;
use num::Float;
use num_traits::cast::FromPrimitive;

fn main() {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();

    // Use Clap to parse command-line arguments
    let matches = Command::new("hnsw2knn")
        .version("0.3.2")
        .about("Extract K Nearest Neighbors (K-NN) from HNSW graph.")
        .arg(
            Arg::new("database_path")
                .short('b')
                .long("hnsw")
                .value_name("DATADIR")
                .help("Directory containing pre-built HNSW database files")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_PATH")
                .help("Output path to write the neighbor list")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("knn")
                .short('n')
                .long("k-nearest-neighbors")
                .value_name("KNN")
                .help("Number of k-nearest-neighbors to extract")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("32"),
        )
        .get_matches();

    // Extract command-line arguments
    let db_path = matches
        .get_one::<String>("database_path")
        .unwrap()
        .to_string();
    let out_path = matches.get_one::<String>("output").unwrap().to_string();
    let knbn = *matches.get_one::<usize>("knn").unwrap();
    // Prepare to reload HNSW from disk
    let database_dirpath = Path::new(&db_path);
    let hnswio_res = reloadhnsw::get_hnswio(database_dirpath);
    if hnswio_res.is_err() {
        panic!("Error: {:?}", hnswio_res.err());
    }
    let mut hnswio = hnswio_res.unwrap();

    // Load sequence dictionary from seqdict.json to retrieve actual sequence names
    let seqdict_path = database_dirpath.join("seqdict.json");
    info!(
        "\nRe-loading sequence dictionary from {}",
        seqdict_path.display()
    );
    let seqdict = SeqDict::reload_json(&seqdict_path);
    let seqdict = match seqdict {
        Ok(seqdict) => seqdict,
        _ => {
            panic!(
                "SeqDict reload from dump file {} failed",
                seqdict_path.display()
            );
        }
    };

    // Load the HNSW graph
    let hnsw_res = hnswio.load_hnsw::<
        <OptDensHashSketch<Kmer32bit,f32> as SeqSketcherT<Kmer32bit>>::Sig,
        DistHamming,
    >();
    if hnsw_res.is_err() {
        panic!("Error: {:?}", hnsw_res.err());
    }
    let hnsw = hnsw_res.unwrap();

    // Choose how many neighbors to extract, the actual neighbors extracted can be small than requested
    // since t also depends on HNSW build maximumlly allowed neighbors
    let kgraph_res = kgraph_from_hnsw_all::<_, _, f32>(&hnsw, knbn);

    // Build the KGraph from HNSW, then save neighbor lists
    match kgraph_res {
        Ok(kgraph) => {
            println!(
                "KGraph successfully created with {} nodes.",
                kgraph.get_nb_nodes()
            );
            // Save the neighbor list to a file, printing actual sequence IDs
            if let Err(e) = save_neighbor_list_to_file(&kgraph, &seqdict, &out_path) {
                eprintln!("Error saving neighbor list: {:?}", e);
            } else {
                println!("Neighbor list saved to {}", out_path);
            }
        }
        Err(e) => {
            eprintln!("Error creating KGraph: {:?}", e);
        }
    }
}

/// Save KGraph neighbor lists to a file, printing *actual sequence IDs* from SeqDict
/// we look up the corresponding file path + FASTA ID from `seqdict`.
fn save_neighbor_list_to_file<F>(
    kgraph: &KGraph<F>,
    seqdict: &SeqDict,
    output_file: &str,
) -> std::io::Result<()>
where
    F: FromPrimitive + Float + std::fmt::UpperExp + Sync + Send + std::iter::Sum,
{
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);

    // Iterate over each node index in the KGraph
    for node_idx in 0..kgraph.get_nb_nodes() {
        // Get the data ID for the current node (this is used to index seqdict)
        if let Some(node_data_id) = kgraph.get_data_id_from_idx(node_idx) {
            // The HNSW data ID is used as an index into seqdict.0[..]
            let node_item = &seqdict.0[*node_data_id];
            let node_path = node_item.get_id().get_path();

            // Write the node's actual sequence ID
            // Format: path|fasta_id:
            write!(writer, "{}:", node_path)?;

            // Get the list of outgoing edges (neighbors)
            let edges = kgraph.get_out_edges_by_idx(node_idx);
            // Iterate over each neighbor
            for edge in edges {
                // Get the neighbor's index
                let neighbor_idx = edge.node;
                // Convert neighbor_idx -> neighbor_data_id -> actual sequence name
                if let Some(neighbor_data_id) = kgraph.get_data_id_from_idx(neighbor_idx) {
                    let neighbor_item = &seqdict.0[*neighbor_data_id];
                    let neighbor_path = neighbor_item.get_id().get_path();

                    // Write neighbor's path|fasta_id plus the edge weight/distance
                    write!(
                        writer,
                        "\t{}:{:.6}",
                        neighbor_path,
                        edge.weight.to_f64().unwrap()
                    )?;
                }
            }
            // Newline after each node's neighbor list
            writeln!(writer)?;
        }
    }

    writer.flush()?;
    Ok(())
}
