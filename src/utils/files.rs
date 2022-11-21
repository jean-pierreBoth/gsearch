//! This file contains directory exploration and fasta file selection

use std::io;
use std::io::{BufReader, BufWriter };

use std::fs::OpenOptions;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};

use rayon::prelude::*;

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
pub fn is_fasta_dna_file(pathb : &PathBuf) -> bool {
    let filename = pathb.to_str().unwrap();
    if filename.ends_with("fna.gz")|| filename.ends_with("fa.gz") || 
                filename.ends_with("fasta.gz") || filename.ends_with("fna") || filename.ends_with("fa") || filename.ends_with("fasta") {
        return true;
    }
    else { 
        log::debug!("found non dna file : {:?}", filename);
        return false;
    }
}  // end of is_fasta_file


/// returns true if file is a fasta file preotein (possibly gzipped) suffixed by .faa
pub fn is_fasta_aa_file(pathb : &PathBuf) -> bool {
    let filename = pathb.to_str().unwrap();
    if filename.ends_with("faa.gz")|| filename.ends_with("faa") {
        return true;
    }
    else { 
        return false;
    }
}  // end of is_fasta_aa_file




/// scan directory recursively, executing function file_task on each file.
/// adapted from from crate fd_find
/// A partially parallel version in process_files_group below
pub fn process_dir(state : &mut ProcessingState, datatype: &DataType, dir: &Path, filter_params : &FilterParams, 
                file_task: &(dyn Fn(&PathBuf, &FilterParams) -> Vec<IdSeq> + Sync), 
                sender : &crossbeam_channel::Sender::<Vec<IdSeq>>) -> io::Result<usize> {
    //
    // we checked that we have a directory
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let pathb = entry.path();
        if pathb.is_dir() {
            let _ =  process_dir(state, &datatype, &pathb, filter_params, file_task, sender)?;
        } else {
            // check if entry is a fasta.gz file or a .faa file
            let mut to_sketch = match datatype {
                DataType::DNA => {
                    if is_fasta_dna_file(&pathb)  {
                        file_task(&pathb, filter_params)
                    }
                    else {
                        log::warn!("process_dir found a non dna file {:?}", entry.file_name());
                        Vec::<IdSeq>::new()
                    }
                },
                DataType::AA => {
                    if is_fasta_aa_file(&pathb) {
                        file_task(&pathb, filter_params)
                    }
                    else {
                        log::warn!("process_dir found a non AA file {:?}", entry.file_name());
                        Vec::<IdSeq>::new()
                    }
                },
            };
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
        } // end of check on datatype
    } // end of for 
    //
    drop(sender);
    //
    Ok(state.nb_seq)
}  // end of process_dirs





/// This function is called by process dirs.
/// aargument entries is a slice of Entry that have been checked to be files (and not directory!)
/// process a block of entry in parallel. Block of entries should be of moderate size. (5 or 10?) depending on the number of threads of the Cpus 
pub fn process_files_group(datatype: &DataType, filter_params : &FilterParams, 
    pathb : &[PathBuf], file_task: &(dyn Fn(&PathBuf, &FilterParams) -> Vec<IdSeq> + Sync)) -> Vec<Vec<IdSeq>> 
{
    //
    let process_file = | pathb : &PathBuf| -> Vec::<IdSeq> {
        let to_sketch = match datatype {
            DataType::DNA => {
                if is_fasta_dna_file(pathb)  {
                    file_task(pathb, filter_params)
                }
                else {
                    log::warn!("process_files_group found a non dna file {:?}", pathb.as_path());
                    Vec::<IdSeq>::new()
                }
            },
            DataType::AA => {
                if is_fasta_aa_file(pathb) {
                    file_task(pathb, filter_params)
                }
                else {
                    log::warn!("process_files_group found a non AA file {:?}", pathb.as_path());
                    Vec::<IdSeq>::new()
                }
            },
        };  
        to_sketch      
        };
    //
    let seqseq = pathb.into_par_iter().map(|file| process_file(file)).collect();
    //
    seqseq
}  // end of process_files_group




pub fn process_dir_parallel(state : &mut ProcessingState, datatype: &DataType, dir: &Path, filter_params : &FilterParams, 
        block_size : usize,
        file_task: &(dyn Fn(&PathBuf, &FilterParams) -> Vec<IdSeq> + Sync), 
        sender : &crossbeam_channel::Sender::<Vec<IdSeq>>) -> io::Result<usize> {
    //
    // We will process files in parallel by blocks of size block_size to control memory
    // and to balance cpu and memory with threads doing hnsw work
    let mut path_block = Vec::<PathBuf>::with_capacity(block_size);
    //
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let pathb = entry.path();
        if pathb.is_dir() {
            let _ =  process_dir(state, &datatype, &pathb, filter_params, file_task, sender)?;
        } 
        else {
            // we have path of a file
            if path_block.len() < block_size {
                path_block.push(pathb.clone());
            }
            else {
                std::panic!("process_dir_parallel buffer overflow , should not occur")
            }
            // if buffer full, process send and empty buffer
            let seqs = process_files_group(datatype,  filter_params, &path_block, file_task);
            for mut seqfile in seqs {
                for i in 0..seqfile.len() {
                    seqfile[i].rank = state.nb_seq;
                    state.nb_seq += 1;
                }                
                sender.send(seqfile).unwrap();
                state.nb_file += 1;
            }
            path_block.clear();
        }
    } // end of for on entry
    // if there is a residue (nbfile not multiple of  BLOCK_SIZE) in path_block we must treat the residue
    if path_block.len() > 0 {
        let seqs = process_files_group(datatype,  filter_params, &path_block, file_task);
        for mut seqfile in seqs {
            for i in 0..seqfile.len() {
                seqfile[i].rank = state.nb_seq;
                state.nb_seq += 1;
            }                
            sender.send(seqfile).unwrap();
            state.nb_file += 1;
        }
        path_block.clear();
    }
    //
    drop(sender);
    //    
    Ok(state.nb_seq)
}  // end of process_dir_parallel