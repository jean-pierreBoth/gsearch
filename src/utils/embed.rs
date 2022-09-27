//! module doing some interface to annembed
//! Provides some statistics on hnw graph reloaded to do archaea identification/  
//! Provide also an embedding if asked for

use annembed::fromhnsw::kgraph_from_hnsw_all;
use annembed::prelude::*;
use csv::Writer;
use hnsw_rs::prelude::*;
use log;
//use ndarray::{s, Array1, Array3};
//use std::fs::OpenOptions;
//use std::io::prelude::*;
//use std::io::BufReader;
//use std::path::PathBuf;

pub fn get_graph_stats_embed<T, D>(hnsw: &Hnsw<T, D>, embed: bool) -> Result<(), ()>
where
    T: Clone + Send + Sync,
    D: Distance<T> + Send + Sync,
{
    let knbn = 8;
    log::info!("calling kgraph_from_hnsw_all, embed = {}", embed);

    let kgraph_res = kgraph_from_hnsw_all::<T, D, f32>(hnsw, knbn);
    if let Ok(kgraph) = kgraph_res {
        // we are just interested in quantile statistics on first distance to neighbours.
        log::info!(" computing graph statistics");
        let _kgraph_stats = kgraph.get_kraph_stats();
        if embed {
            log::info!(" going to embedding");
            let mut embed_params = EmbedderParams::default();
            embed_params.nb_grad_batch = 15;
            embed_params.scale_rho = 0.5;
            embed_params.beta = 1.;
            embed_params.grad_step = 3.;
            embed_params.nb_sampling_by_edge = 10;
            embed_params.dmap_init = true;
            let mut embedder = Embedder::new(&kgraph, embed_params);
            let embed_res = embedder.embed();

            if embed_res.is_err() {
                log::error!("embedding failed");
                std::process::exit(1);
            }
            //dump embeddings into csv file
            log::info!("dumping in csv file");
            let mut csv_w = Writer::from_path("database_embedded.csv").unwrap();
            // we can use get_embedded_reindexed as we indexed DataId contiguously in hnsw!
            let _res = write_csv_array2(&mut csv_w, &embedder.get_embedded_reindexed());
            csv_w.flush().unwrap();

            log::info!(" embedding finished");
        } // end of embedding
    } else {
        log::error!("could not get graph stats from hnsw");
        return Err(());
    }
    Ok(())
} // end of get_graph_stats
