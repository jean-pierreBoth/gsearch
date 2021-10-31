//! module doing some interface to annembed
//! Provides some statistics on hnw graph reloaded to do archaea identification/  
//! Provide also an embedding if asked for

use hnsw_rs::prelude::*;

use annembed::prelude::*;

pub fn get_graph_stats_embed<T, D>(hnsw: &Hnsw<T,D>, embed : bool) -> Result<(),()>
    where T : Clone+Send+Sync,  
          D : Distance<T> + Send + Sync {
    let knbn = 8;
    log::info!("calling kgraph_from_hnsw_all");

    let kgraph_res = kgraph_from_hnsw_all::<T, D, f32>(hnsw, knbn);
    if let Ok(kgraph) = kgraph_res {
        // we are just interested in quantile statistics on first distance to neighbours.
        let _kgraph_stats = kgraph.get_kraph_stats();
        if embed {
            let mut embed_params = EmbedderParams::new();
            embed_params.nb_grad_batch = 15;
            embed_params.scale_rho = 0.5;
            embed_params.beta = 1.;
            embed_params.grad_step = 3.;
            embed_params.nb_sampling_by_edge = 10;
            embed_params.dmap_init = true;
            let mut embedder = Embedder::new(&kgraph, embed_params);
            let _embed_res = embedder.embed();
        } // end of embedding

    }
    else {
        log::error!("could not get graph stats from hnsw");
        return Err(())
    }
    Ok(())
}  // end of get_graph_stats


