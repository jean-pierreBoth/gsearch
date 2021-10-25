//! structures related to processing parameters


use std::fs::OpenOptions;
use std::path::{Path};
use std::io::{BufReader, BufWriter };

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};


use kmerutils::sketching::seqsketchjaccard::SeqSketcher;

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


/// Gathers parameters used for hnsw, sketching and choice of sequence/blocked genome processing.
/// To be dumped to ease request processing.
#[derive(Copy,Clone,Serialize,Deserialize)]
pub struct ProcessingParams {
    /// hnsw parameters
    hnsw: HnswParams,
    ///
    sketch: SeqSketcher,
    /// do we process sequence by sequence? true if we process whole genome in one block, false if we process sequence by sequence 
    block_flag : bool,
}  // end of WholeParams



impl ProcessingParams {

    pub fn new(hnsw : HnswParams, sketch :  SeqSketcher, block_flag : bool) -> Self {
        ProcessingParams{hnsw, sketch, block_flag}
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



    /// reload from a json dump
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

