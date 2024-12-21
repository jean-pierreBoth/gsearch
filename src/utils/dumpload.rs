//! a small module dedicated to dump reload parameters, seqdict ...
//! common to dna and aa processing
//!

use serde::{de::DeserializeOwned, Serialize};
use std::path::{Path, PathBuf};

use std::fmt::Debug;

use hnsw_rs::prelude::*;

use crate::utils::{idsketch::*, parameters::*};

// This function dumps hnsw , seqdict and processing params in the same directory given by dump_path_ref
pub fn dumpall<Sig>(
    dump_path_ref: &PathBuf,
    hnsw: &Hnsw<Sig, DistHamming>,
    seqdict: &SeqDict,
    processing_params: &ProcessingParams,
) -> anyhow::Result<(), String>
where
    Sig: Clone + Copy + Send + Sync + Serialize + DeserializeOwned + Debug,
    DistHamming: Distance<Sig>,
{
    if hnsw.get_nb_point() > 0 {
        let mut hnsw_dump = dump_path_ref.to_path_buf().clone();
        hnsw_dump.push("hnswdump");
        let hnswdumpname = String::from(hnsw_dump.to_str().unwrap());
        log::info!("going to dump hnsw with prefix : {:?}", hnswdumpname);
        let dir = PathBuf::from(".");
        let resdump = hnsw.file_dump(&dir, &hnswdumpname);
        match resdump {
            Err(msg) => {
                println!("dump failed error msg : {}", msg);
            }
            _ => {
                println!("dump of hnsw ended");
            }
        };
        // dump some info on layer structure
        hnsw.dump_layer_info();
        // dumping dictionary
        let mut seq_pb = dump_path_ref.clone();
        seq_pb.push("seqdict.json");
        let seqdict_name = String::from(seq_pb.to_str().unwrap());
        let resdump = seqdict.dump(seqdict_name);
        match resdump {
            Err(msg) => {
                println!("seqdict dump failed error msg : {}", msg);
            }
            _ => {
                println!("dump of seqdict ended OK");
            }
        };
    } else {
        log::info!("no dumping hnsw, no data points");
    }
    // and finally dump processing parameters in file name "parameters.json"
    let _ = processing_params.dump_json(dump_path_ref);
    //
    Ok(())
} // end of dumpall

pub fn reload_seqdict(dump_path_ref: &PathBuf) -> SeqDict {
    // must reload seqdict
    let mut filepath = dump_path_ref.clone();
    filepath.push("seqdict.json");
    let res_reload = SeqDict::reload_json(&filepath);
    if res_reload.is_err() {
        let cwd = std::env::current_dir();
        if cwd.is_ok() {
            log::info!("current directory : {:?}", cwd.unwrap());
        }
        log::error!("cannot reload SeqDict (file 'seq.json' from current directory");
        std::process::exit(1);
    } else {
        
        res_reload.unwrap()
    }
} // end of reload_seqdict

/// reload json dumped SeqDict and re-dump a summary a csv file containing filename.
// useful for clustering (coreset)
pub fn seqdict_jsontocsv(json_path_ref: &PathBuf) -> usize {
    println!("in seqdict_jsontocsv");
    let seqdict = reload_seqdict(json_path_ref);
    //
    let mut seqdict_csv_path = json_path_ref.clone();
    seqdict_csv_path.push(Path::new("seqdict.csv"));
    log::info!(
        "seqdict_jsontocsv rdump of seqdict in csv file : {:?}",
        seqdict_csv_path
    );
    //
    let rec_count = seqdict.dump_csv(&seqdict_csv_path);
    log::info!("seqdict_jsontocsv dumped nb record : {}", rec_count);
    //
    rec_count
}

// retrieve or allocate a SeqDict depending on use case
pub fn get_seqdict(
    dump_path_ref: &PathBuf,
    other_params: &ComputingParams,
) -> anyhow::Result<SeqDict> {
    //
    let seqdict = if other_params.get_adding_mode() {
        reload_seqdict(dump_path_ref)
    } else {
        SeqDict::new(100000)
    };
    //
    Ok(seqdict)
} // end of get_seqdict
