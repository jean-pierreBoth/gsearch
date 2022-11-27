//! structures related to processing parameters


use std::fs::OpenOptions;
use std::path::{Path};
use std::io::{BufReader, BufWriter };

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};



/// a structure to filter files or sequences we treat
pub struct FilterParams {
    /// minimum sequence size
    pub min_seq_size : usize,
} // end of struct FilterParams


impl FilterParams {
    pub fn new(min_seq_size : usize) -> Self {
        FilterParams{min_seq_size}
    } // end of new

    /// returns true if we filter (garbage the sequence)
    pub fn filter(&self, seq : &[u8]) -> bool {
        if seq.len() < self.min_seq_size {
            true
        }
        else {
            false
        }
    }
}  // end of FilterParams

//===========================================================

#[derive(Copy,Clone,Serialize,Deserialize)]
pub struct HnswParams {
    /// expected number of sequences to store
    pub capacity : usize,
    /// ef used in construction
    ef : usize,
    max_nb_conn : u8,
}

impl HnswParams {
    pub fn new(capacity : usize, ef : usize, max_nb_conn : u8) -> Self {
        HnswParams{capacity, ef, max_nb_conn}
    }
    //
    pub fn get_ef(&self) -> usize {
        return self.ef
    }

    pub fn get_max_nb_connection(&self) -> u8 {
        self.max_nb_conn
    }
}  // end of impl block HnswParams


//======================================================================================

#[derive(Clone, Copy, Serialize,Deserialize)]
pub struct AnnParameters {
    /// true to get statistics on neighbours
    ask_stats : bool,
    /// do an embedding or not
    embed : bool,
} // end of 


impl Default for AnnParameters {
    fn default() -> Self {
        AnnParameters{ask_stats:false, embed:false}
    }
} // end of default for AnnParameters


impl AnnParameters {
    // default initiaization is false for embed
    pub fn new(ask_stats : bool, embed : bool) -> Self {
        AnnParameters{ask_stats, embed}
    }

    /// returns if stats on distance on nearest neighbours resulting from hnsw information were asked for
    /// In this case some quantiles on nearest neighbours are dumped.
    pub fn ask_stats(&self) -> bool {
        self.ask_stats
    }

    pub fn embed(&self) -> bool {
        self.embed
    }
} // end of impl AnnParameters


//=========================================================================================

// This is redundant with struct Sketcher for DNA case and RNA case, but it makes
// possible the factorization of all parameters

#[derive(Copy,Clone,Serialize,Deserialize)]
pub struct SeqSketcherParams {
    kmer_size : usize,
    sketch_size : usize  
}


impl SeqSketcherParams {
    /// 
    pub fn new(kmer_size: usize, sketch_size : usize) -> Self {
        SeqSketcherParams{kmer_size, sketch_size}
    }

    /// returns kmer size
    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /// return sketch size
    pub fn get_sketch_size(&self) -> usize {
        self.sketch_size
    }  

}  // end of SeqSketcherParams

//==========================================================================================

/// Gathers parameters used for hnsw, sketching and choice of sequence/blocked genome processing.
/// To be dumped to ease request processing.
#[derive(Copy,Clone,Serialize,Deserialize)]
pub struct ProcessingParams {
    /// hnsw parameters
    hnsw: HnswParams,
    /// sketching prameters
    sketch: SeqSketcherParams,
    /// do we process sequence by sequence? true if we process whole genome in one block, false if we process sequence by sequence 
    block_flag : bool,
}  // end of WholeParams



impl ProcessingParams {

    pub fn new(hnsw : HnswParams, sketch :  SeqSketcherParams, block_flag : bool) -> Self {
        ProcessingParams{hnsw, sketch, block_flag}
    }

    /// get parameters for hnsw construction
    pub fn get_hnsw_params(&self) -> &HnswParams {
        &self.hnsw
    }

    /// get parameters used for sketching sequences 
    pub fn get_sketching_params(&self) -> &SeqSketcherParams {
        &self.sketch
    }

    /// get block flag
    pub fn get_block_flag(&self) -> bool {
        self.block_flag
    }

    /// return kmer size used in sketching.
    pub fn get_kmer_size(&self) ->  usize {
        self.sketch.get_kmer_size()
    }


    pub fn dump_json(&self, dirpath: &Path) ->  Result<(), String> {
        //
        let filepath = dirpath.join("parameters.json");
        //
        log::info!("dumping ProcessingParams in json file : {:?}", filepath);
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("ProcessingParams dump : dump could not open file {:?}", filepath.as_os_str());
            println!("ProcessingParams dump: could not open file {:?}", filepath.as_os_str());
            return Err("ProcessingParams dump failed".to_string());
        }
        // 
        let mut writer = BufWriter::new(fileres.unwrap());
        let _ = to_writer(&mut writer, &self).unwrap();
        //
        Ok(())
    } // end of dump_json



    /// reload from a json dump. Used in request module to ensure coherence with database constitution
    pub fn reload_json(dirpath : &Path) -> Result<Self, String> {
        log::info!("in reload_json");
        //
        let filepath = dirpath.join("parameters.json");
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!("ProcessingParams reload_json : reload could not open file {:?}", filepath.as_os_str());
            println!("ProcessingParams reload_json: could not open file {:?}", filepath.as_os_str());
            return Err("ProcessingParams reload_json could not open file".to_string());            
        }
        //
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        let processing_state:Self = serde_json::from_reader(reader).unwrap();
        //
        log::info!("ProcessingParameters reload, blocked : {}", processing_state.block_flag);     
        //
        Ok(processing_state)
    } // end of reload_json


} // end of impl ProcessingParams


//=====================================================================================

/// Some others parameters not necessary for reload
pub struct ComputingParams {
    /// set to true when parsing fasta files in done in //
    parallel_io : bool,
    /// set to true when we increase a hnsw database
    adding_mode : bool,
}


impl Default for ComputingParams {
    fn default() -> Self {
        ComputingParams{parallel_io: false , adding_mode: false}
    }
}


impl ComputingParams {
    pub fn new(parallel_io : bool, adding_mode : bool) -> Self {
        ComputingParams{parallel_io, adding_mode}
    }

    pub fn get_parallel_io(&self) -> bool {
        self.parallel_io
    }

    pub fn get_adding_mode(&self) -> bool {
        self.adding_mode
    }
} // end of ComputingParams