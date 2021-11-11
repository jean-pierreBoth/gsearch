//! This file contains directory exploration and fasta file selection

use std::io;
use std::io::{BufReader, BufWriter };
use std::fs::{self, DirEntry};

use std::fs::OpenOptions;
use std::path::{Path};

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};



use super::idsketch::{IdSeq};
use super::parameters::*;


/// To keep track of processed file and sequence processed
#[derive(Serialize,Deserialize, Clone, Copy)]
pub struct ProcessingState {
    /// nb sequences processed
    pub nb_seq : usize,
    /// nb file processed
    pub nb_file : usize,
    /// elapsed time in sec
    pub elapsed_t : f32,
} 


impl ProcessingState {
    pub fn new() -> Self {
        ProcessingState{nb_seq : 0, nb_file : 0, elapsed_t : 0.}
    } // end of new
    
    /// get elapsed time
    pub fn get_elapsed_t(&self) -> f32 {
        self.elapsed_t
    }
    /// serialized dump
    pub fn dump_json(&self, dirpath : &Path) -> Result<(), String> {
        //
        let filepath = dirpath.join("processing_state.json");
        //
        log::info!("dumping ProcessingState in json file : {:?}", filepath);
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("ProcessingState dump : dump could not open file {:?}", filepath.as_os_str());
            println!("ProcessingState dump: could not open file {:?}", filepath.as_os_str());
            return Err("ProcessingState dump failed".to_string());
        }
        // 
        let mut writer = BufWriter::new(fileres.unwrap());
        let _ = to_writer(&mut writer, &self).unwrap();
        //
        Ok(())
    } // end of dump



    /// reload from a json dump
    pub fn reload_json(dirpath : &Path) -> Result<Self, String> {
        log::info!("in reload_json");
        //
        let filepath = dirpath.join("processing_state.json");
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!("ProcessingState reload_json : reload could not open file {:?}", filepath.as_os_str());
            println!("ProcessingState reload_json: could not open file {:?}", filepath.as_os_str());
            return Err("ProcessingState reload_json could not open file".to_string());            
        }
        //
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        let processing_state:Self = serde_json::from_reader(reader).unwrap();
        //
        log::info!("ProcessingState reload, nb sequences : {}", processing_state.nb_seq);     
        //
        Ok(processing_state)
    } // end of reload_json
} // end of ProcessingState



///
/// enum to store if we are doing DNA or AA sequence processing mode
pub enum DataType {
    DNA,
    AA,
}

//==================================================================================


// returns true if file is a fasta file (possibly gzipped)
// filename are of type GCA[GCF]_000091165.1_genomic.fna.gz
pub fn is_fasta_dna_file(file : &DirEntry) -> bool {
    let filename = file.file_name().into_string().unwrap();
    if filename.ends_with("fna.gz")|| filename.ends_with("fa.gz") || filename.ends_with("fasta.gz") {
        return true;
    }
    else { 
        return false;
    }
}  // end of is_fasta_file


/// returns true if file is a fasta file RNA (possibly gzipped) suffixed by .faa
pub fn is_fasta_rna_file(file : &DirEntry) -> bool {
    let filename = file.file_name().into_string().unwrap();
    if filename.ends_with("faa.gz")|| filename.ends_with("faa") {
        return true;
    }
    else { 
        return false;
    }
}  // end of is_fasta_aa_file











/// scan directory recursively, executing function file_task on each file.
/// adapted from from crate fd_find
pub fn process_dir(state : &mut ProcessingState, dir: &Path, filter_params : &FilterParams, 
                file_task: &dyn Fn(&DirEntry, &FilterParams) -> Vec<IdSeq>, 
                sender : &crossbeam_channel::Sender::<Vec<IdSeq>>) -> io::Result<usize> {
    //
    // we checked that we have a directory
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            let _ =  process_dir(state, &path, filter_params, file_task, sender)?;
        } else {
            // check if entry is a fasta.gz file or a .faa file
            // TODO should check that there is no mix of files?
            if is_fasta_dna_file(&entry) || is_fasta_rna_file(&entry) {
                let mut to_sketch = file_task(&entry, filter_params);
                // put a rank id in sequences, now we have full information of where do the sequence come from
                for i in 0..to_sketch.len() {
                    to_sketch[i].rank = state.nb_seq;
                    state.nb_seq += 1;
                }
                state.nb_file += 1;
                if log::log_enabled!(log::Level::Info) && state.nb_file % 500 == 0 {
                    log::info!("nb file processed : {}, nb sequences processed : {}", state.nb_file, state.nb_seq);
                }
                if state.nb_file % 500 == 0 {
                    println!("nb file processed : {}, nb sequences processed : {}", state.nb_file, state.nb_seq);
                }
                // we must send to_sketch into channel to upper thread
                sender.send(to_sketch).unwrap();
            }
        }
    }
    //
    drop(sender);
    //
    Ok(state.nb_seq)
}  // end of visit_dirs

