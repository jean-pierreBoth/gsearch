//! module doing some interface to annembed
//! Provides some statistics on hnw graph reloaded to do genome identification/  
//! Provide also an embedding if asked for
//! 

use csv::Writer;

// our use
use hnsw_rs::prelude::*;
use annembed::fromhnsw::{hubness::Hubness,kgraph_from_hnsw_all};
use annembed::fromhnsw::kgraph::KGraph;
use annembed::prelude::*;
use cpu_time::ProcessTime;
use std::time::SystemTime;
use std::time::Duration;

pub fn get_graph_stats_embed<T, D>(hnsw: &Hnsw<T, D>, embed: bool, embed_params_opt : Option<EmbedderParams>) -> Result<(), ()>
where
    T: Clone + Send + Sync,
    D: Distance<T> + Send + Sync,
{
    let knbn = 15;
    log::info!("calling kgraph_from_hnsw_all, embed = {}", embed);

    let kgraph_res = kgraph_from_hnsw_all::<T, D, f32>(hnsw, knbn);
    if let Ok(kgraph) = kgraph_res {
        // we are just interested in quantile statistics on first distance to neighbours.
        log::info!("\n\n computing graph statistics");
        let _kgraph_stats = kgraph.get_kraph_stats();
        
        // Local Intrinsic dimension and hubness
        // Get some statistics on induced graph. This is not related to the embedding process
        let knbn_new = 25;
        let kgraph_new : KGraph<f32>;
        kgraph_new = kgraph_from_hnsw_all(&hnsw, knbn_new).unwrap();
        log::info!("minimum number of neighbours {}", kgraph_new.get_max_nbng());
        log::info!("dimension estimation...");
        let sampling_size = 10000;
        let cpu_start = ProcessTime::now();
        let sys_now = SystemTime::now();
        let dim_stat = kgraph_new.estimate_intrinsic_dim(sampling_size);
        let cpu_time: Duration = cpu_start.elapsed();
        println!("\n dimension estimation sys time(ms) : {:.3e},  cpu time(ms) {:.3e}\n", sys_now.elapsed().unwrap().as_millis(), cpu_time.as_millis());
        if dim_stat.is_ok() {
            let dim_stat = dim_stat.unwrap();
            log::info!("\n dimension estimation with nbpoints : {}, dim : {:.3e}, sigma = {:.3e} \n", 
                sampling_size, dim_stat.0, dim_stat.1);
            println!(" dimension estimation with nbpoints : {}, dim : {:.3e}, sigma = {:.3e}", 
                sampling_size, dim_stat.0, dim_stat.1); 
        }
        
        log::info!("\n\n hubness summary");
        let hubness = Hubness::new(&kgraph);
        let s3_hubness = hubness.get_standard3m();
        log::info!("\n graph hubness estimation : {:.3e}", s3_hubness);
        println!("\n graph hubness estimation : {:.3e} \n", s3_hubness);
        let _ = hubness.get_hubness_histogram();
        //
        if embed {
            log::info!(" going to embedding");
            let mut embed_params = match embed_params_opt {
                Some(embed_params) => { embed_params },
                None               => {
                                        let mut embed_params = EmbedderParams::default();
                                        embed_params.nb_grad_batch = 15;
                                        embed_params.scale_rho = 0.75;
                                        embed_params.beta = 1.;
                                        embed_params.grad_step = 3.;
                                        embed_params.nb_sampling_by_edge = 10;
                                        embed_params.dmap_init = true;
                                        embed_params
                }
            };
            //
            if hnsw.get_point_indexation().get_layer_nb_point(1) > 30000 {
                log::info!("doing hierarchical embedding from layer 1");
                embed_params.set_hierarchy_layer(1);
            }
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
            // try to get a quality indicator
            log::info!("calling embedder.get_quality_estimate_from_edge_length");
            let _quality = embedder.get_quality_estimate_from_edge_length(200);
        } // end of embedding
    } else {
        log::error!("could not get graph stats from hnsw");
        return Err(());
    }
    Ok(())
} // end of get_graph_stats
