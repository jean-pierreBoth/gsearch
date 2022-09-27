//! Reload the hnsw structure with json parameters



use std::path::{Path};
use std::fs::{OpenOptions};
use std::io::{BufReader};
use serde::{Serialize, de::DeserializeOwned};

use std::time::{SystemTime};

use hnsw_rs::prelude::*;
use hnsw_rs::hnswio::{load_description, load_hnsw};

//mod files;
use crate::utils::AnnParameters;



/* 
    reload hnsw from dump directory
    We know filename : hnswdump.hnsw.data and hnswdump.hnsw.graph
    The allow(unused_variables) is here to avoid a warning on ann_params when compiling without feature f_annembed
 */
pub fn reload_hnsw<T>(dump_dirpath : &Path, _ann_params: &AnnParameters) -> Option<Hnsw<T, DistHamming>>  
            where T : 'static + Clone + Send + Sync + Serialize + DeserializeOwned ,
                DistHamming : Distance<T>  {
    // just concat dirpath to filenames and get pathbuf
    let graph_path = dump_dirpath.join("hnswdump.hnsw.graph");
    log::info!("reload_hnsw, loading graph from {}",graph_path.to_str().unwrap());
    let graphfile = OpenOptions::new().read(true).open(&graph_path);
    if graphfile.is_err() {
        println!("test_dump_reload: could not open file {:?}", graph_path.as_os_str());
        return None;
    }
    let graphfile = graphfile.unwrap();
    let mut graphfile = BufReader::with_capacity(50_000_000, graphfile);
    //
    let data_path = dump_dirpath.join("hnswdump.hnsw.data");
    log::info!("reload_hnsw, loading data from {}",data_path.to_str().unwrap());
    let datafile = OpenOptions::new().read(true).open(&data_path);
    if datafile.is_err() {
        println!("test_dump_reload: could not open file {:?}", data_path.as_os_str());
        return None;
    }
    let datafile = datafile.unwrap();
    let mut datafile = BufReader::with_capacity(50_000_000,datafile);
    //
    let start_t = SystemTime::now();
    let hnsw_description = load_description(&mut graphfile).unwrap();
    let hnsw : Hnsw<T, DistHamming>= load_hnsw(&mut graphfile, &hnsw_description, &mut datafile).unwrap();
    let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
    if log::log_enabled!(log::Level::Info) {
        log::info!("reload_hnsw : elapsed system time(s) {}", elapsed_t);
    }
    else {
        println!("reload_hnsw : elapsed system time(s) {}", elapsed_t);
    }
    // feature enabled (or not) in Cargo.toml, requires the crate annembed
    #[cfg(any(feature="annembed_openblas-system", feature="annembed_openblas-static" , feature="annembed_intel-mkl"))]
    if _ann_params.ask_stats() {
        let _ = super::embed::get_graph_stats_embed(&hnsw, true);
    }
    //
    return Some(hnsw);
    //  
} // end of reload_hnsw

