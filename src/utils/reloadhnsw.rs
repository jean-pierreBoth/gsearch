//! Reload the hnsw structure with json parameters



use std::path::Path;
use std::fs::OpenOptions;
use std::io::BufReader;


use hnsw_rs::hnswio::*;

#[allow(unused)]
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


///    reload hnsw from dump directory
///    We know filename : hnswdump.hnsw.data and hnswdump.hnsw.graph
pub fn get_hnswio(dump_dirpath : &Path) -> Result<HnswIo, String> {
    //
    let dumpname = String::from("hnswdump.hnsw");
    let pathb = std::path::PathBuf::from(dump_dirpath);
    let mut reloader = HnswIo::new(pathb, dumpname);
    let options = ReloadOptions::default();
    reloader.set_options(options);
    //
    return Ok(reloader);
    //  
} // end of reload_hnsw

