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


// we must know data type used in dump!
pub fn get_hnsw_type(dump_dirpath : &Path) -> anyhow::Result<String>  {
    let graph_path = dump_dirpath.join("hnswdump.hnsw.graph");
    log::info!("reload_hnsw, loading graph from {}",graph_path.to_str().unwrap());
    let graphfile = OpenOptions::new().read(true).open(&graph_path);
    if graphfile.is_err() {
        println!("reload_hnsw: could not open file {:?}", graph_path.as_os_str());
        return Err(anyhow::anyhow!("reload_hnsw: could not open file"));
    }
    let graphfile = graphfile.unwrap();
    let mut graphfile = BufReader::with_capacity(50_000, graphfile);
    let hnsw_description = load_description(&mut graphfile);
    if hnsw_description.is_err() {
        log::error!("reload_hnsw: could not reload description");
        return Err(anyhow::anyhow!("reload_hnsw: could not reload description"));
    }
    let hnsw_description = hnsw_description.unwrap();
    //
    return Ok(hnsw_description.t_name.clone());
}  // end of get_hnsw_type


/* 
    reload hnsw from dump directory
    We know filename : hnswdump.hnsw.data and hnswdump.hnsw.graph
 */
pub fn reload_hnsw<T>(dump_dirpath : &Path, _ann_params: &AnnParameters) -> Result<Hnsw<T, DistHamming>, String>  
            where T : 'static + Clone + Send + Sync + Serialize + DeserializeOwned ,
                DistHamming : Distance<T>  {
    // just concat dirpath to filenames and get pathbuf
    let graph_path = dump_dirpath.join("hnswdump.hnsw.graph");
    log::info!("reload_hnsw, loading graph from {}",graph_path.to_str().unwrap());
    let graphfile = OpenOptions::new().read(true).open(&graph_path);
    if graphfile.is_err() {
        println!("reload_hnsw: could not open file {:?}", graph_path.as_os_str());
        return Err(String::from("reload_hnsw: could not open file"));
    }
    let graphfile = graphfile.unwrap();
    let mut graphfile = BufReader::with_capacity(50_000_000, graphfile);
    //
    let data_path = dump_dirpath.join("hnswdump.hnsw.data");
    log::info!("reload_hnsw, loading data from {}",data_path.to_str().unwrap());
    let datafile = OpenOptions::new().read(true).open(&data_path);
    if datafile.is_err() {
        println!("test_dump_reload: could not open file {:?}", data_path.as_os_str());
        return Err(String::from("reload_hnsw: could not open file"));
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
    if _ann_params.ask_stats() || _ann_params.embed() {
        log::info!("calling embed::get_graph_stats_embed");
        let _ = super::embed::get_graph_stats_embed(&hnsw, _ann_params.embed());
    }
    //
    return Ok(hnsw);
    //  
} // end of reload_hnsw

