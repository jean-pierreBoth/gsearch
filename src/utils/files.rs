//! This file contains directory exploration and fasta file selection

use std::io;
use std::io::{BufReader, BufWriter};

use std::fs::{File, OpenOptions};
use std::io::Read;
use std::path::{Path, PathBuf};

use cpu_time::ThreadTime;
use std::time::SystemTime;

use serde::{Deserialize, Serialize};
use serde_json::to_writer;

use rayon::prelude::*;

use super::idsketch::IdSeq;
use super::parameters::*;

/// To keep track of processed file and sequence processed
#[derive(Serialize, Deserialize, Clone, Copy)]
pub struct ProcessingState {
    /// nb sequences processed
    pub nb_seq: usize,
    /// nb file processed
    pub nb_file: usize,
    /// elapsed time in sec
    pub elapsed_t: f32,
}

impl Default for ProcessingState {
    fn default() -> Self {
        Self::new()
    }
}

impl ProcessingState {
    pub fn new() -> Self {
        ProcessingState {
            nb_seq: 0,
            nb_file: 0,
            elapsed_t: 0.,
        }
    } // end of new

    /// get elapsed time
    pub fn get_elapsed_t(&self) -> f32 {
        self.elapsed_t
    }
    /// serialized dump
    pub fn dump_json(&self, dirpath: &Path) -> Result<(), String> {
        //
        let filepath = dirpath.join("processing_state.json");
        //
        log::info!("dumping ProcessingState in json file : {:?}", filepath);
        //
        let fileres = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&filepath);
        if fileres.is_err() {
            log::error!(
                "ProcessingState dump : dump could not open file {:?}",
                filepath.as_os_str()
            );
            println!(
                "ProcessingState dump: could not open file {:?}",
                filepath.as_os_str()
            );
            return Err("ProcessingState dump failed".to_string());
        }
        //
        let mut writer = BufWriter::new(fileres.unwrap());
        to_writer(&mut writer, &self).unwrap();
        //
        Ok(())
    } // end of dump

    /// reload from a json dump
    pub fn reload_json(dirpath: &Path) -> Result<Self, String> {
        //
        let filepath = dirpath.join("processing_state.json");
        log::info!("reloading processing state from : {:?}", filepath);
        //
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!(
                "ProcessingState reload_json : reload could not open file {:?}",
                filepath.as_os_str()
            );
            println!(
                "ProcessingState reload_json: could not open file {:?}",
                filepath.as_os_str()
            );
            return Err("ProcessingState reload_json could not open file".to_string());
        }
        //
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        let processing_state: Self = serde_json::from_reader(reader).unwrap();
        //
        log::info!(
            "ProcessingState reload, nb sequences : {}",
            processing_state.nb_seq
        );
        //
        Ok(processing_state)
    } // end of reload_json
} // end of ProcessingState

//==================================================================================

// returns true if file is a fasta file (possibly gzipped)
// filename are of type GCA[GCF]_000091165.1_genomic.fna.gz
pub fn is_fasta_dna_file(pathb: &PathBuf) -> bool {
    let filename = pathb.to_str().unwrap();
    if filename.ends_with("fna.gz")
        || filename.ends_with("fa.gz")
        || filename.ends_with("fa.xz")
        || filename.ends_with("fna.xz")
        || filename.ends_with("fasta.xz")
        || filename.ends_with("fa.bz2")
        || filename.ends_with("fna.bz2")
        || filename.ends_with("fasta.bz2")
        || filename.ends_with("fasta.gz")
        || filename.ends_with("fna")
        || filename.ends_with("fa")
        || filename.ends_with("fasta")
    {
        true
    } else {
        log::debug!("found non dna file : {:?}", filename);
        false
    }
} // end of is_fasta_file

/// returns true if file is a fasta file preotein (possibly gzipped) suffixed by .faa
pub fn is_fasta_aa_file(pathb: &PathBuf) -> bool {
    let filename = pathb.to_str().unwrap();
    filename.ends_with("faa.gz")
        || filename.ends_with("faa")
        || filename.ends_with("faa.xz")
        || filename.ends_with("faa.bz2")
} // end of is_fasta_aa_file

/// scan directory recursively, executing function file_task on each file.
/// adapted from from crate fd_find
/// A partially parallel version in process_files_group below
pub fn process_dir(
    state: &mut ProcessingState,
    datatype: &DataType,
    dir: &Path,
    filter_params: &FilterParams,
    file_task: &dyn Fn(&PathBuf, &FilterParams) -> Vec<IdSeq>,
    sender: &crossbeam_channel::Sender<Vec<IdSeq>>,
) -> io::Result<usize> {
    //
    // we checked that we have a directory
    let nb_files_at_entry = state.nb_file;
    let mut nb_sent = 0;
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let pathb = entry.path();
        if pathb.is_dir() {
            nb_sent = process_dir(state, datatype, &pathb, filter_params, file_task, sender)?;
        } else {
            // check if entry is a fasta.gz file or a .faa file
            let mut to_sketch = match datatype {
                DataType::DNA => {
                    if is_fasta_dna_file(&pathb) {
                        Some(file_task(&pathb, filter_params))
                    } else {
                        log::warn!("process_dir found a non dna file {:?}", entry.file_name());
                        None
                    }
                }
                DataType::AA => {
                    if is_fasta_aa_file(&pathb) {
                        Some(file_task(&pathb, filter_params))
                    } else {
                        log::warn!("process_dir found a non AA file {:?}", entry.file_name());
                        None
                    }
                }
            };
            // put a rank id in sequences, now we have full information of where do the sequence come from
            if to_sketch.is_some() {
                let to_sketch_ref = to_sketch.as_mut().unwrap();
                for i in 0..to_sketch_ref.len() {
                    to_sketch_ref[i].rank = state.nb_seq;
                    to_sketch_ref[i].filerank = state.nb_file;
                    state.nb_seq += 1;
                }
                state.nb_file += 1;
                if log::log_enabled!(log::Level::Info) && state.nb_file % 1000 == 0 {
                    log::info!(
                        "nb file processed : {}, nb sequences processed : {}",
                        state.nb_file,
                        state.nb_seq
                    );
                }
                if state.nb_file % 1000 == 0 {
                    println!(
                        "nb file processed : {}, nb sequences processed : {}",
                        state.nb_file, state.nb_seq
                    );
                }
                // we must send to_sketch into channel to upper thread
                sender.send(to_sketch.unwrap()).unwrap();
            }
        } // end of check on datatype
    } // end of for 
    //
    Ok(nb_sent + state.nb_file - nb_files_at_entry)
} // end of process_dirs

/// open, just read whole file and return buffer for further Fasta processing by needletail
fn file_to_buffer(pathb: &PathBuf) -> Vec<u8> {
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let metadata = std::fs::metadata(pathb.clone());
    let f_len: usize;
    match metadata {
        Ok(metadata) => {
            f_len = metadata.len() as usize;
        }
        Err(_) => {
            log::error!(
                "file_to_buffer could not get length of file {} ",
                pathb.to_str().unwrap()
            );
            f_len = 10_000_000;
        }
    }
    let mut readfile: Vec<u8> = Vec::with_capacity(f_len);
    let mut bufread = std::io::BufReader::with_capacity(f_len + 1000, File::open(pathb).unwrap());
    let res = bufread.read_to_end(&mut readfile);
    assert!(res.is_ok());
    let nb_read = res.unwrap();
    log::trace!(
        "nb byte read from : {:?} , {}, length : {}",
        pathb.clone(),
        nb_read,
        f_len
    );
    //
    readfile
} // end of file_to_buffer

/// This function is called by process_dir_parallel_rec
/// aargument entries is a slice of directory Entry that have been checked to be files (and not directory!).  
/// Files have ben previously sequentially read and transformed into u8 slices, then we can execute file_task in parallel
/// slices associated to files. fIn our usage ile_task does the decompressing and parsing of fasta files.  
///
/// **Block of entries should be of size depending on the number of threads of the Cpus and the memory at disposal and the size of decompressed files**.
///
pub(crate) fn process_files_group(
    datatype: &DataType,
    filter_params: &FilterParams,
    pathb: &[PathBuf],
    file_task: &(dyn Fn(&PathBuf, &[u8], &FilterParams) -> Vec<IdSeq> + Sync),
) -> Vec<Vec<IdSeq>> {
    //
    if log::log_enabled!(log::Level::Debug) {
        log::debug!("process_files_group recieved : {} files", pathb.len());
        log::debug!("memory  : {:?}", memory_stats::memory_stats().unwrap());
    }
    //
    let start_t = SystemTime::now();
    let cpu_start = ThreadTime::try_now();
    //
    let process_file = |pathb: &PathBuf, buffer: &[u8]| -> Option<Vec<IdSeq>> {
        if !buffer.is_empty() {
            Some(file_task(pathb, buffer, filter_params))
        } else {
            None
        }
    };
    // get read results for each files
    let files_read: Vec<Vec<u8>> = pathb
        .iter()
        .map(|path| match datatype {
            DataType::DNA => {
                if is_fasta_dna_file(path) {
                    file_to_buffer(path)
                } else {
                    log::warn!(" encountering a not dna file : {:?}", path);
                    Vec::<u8>::new()
                }
            }
            DataType::AA => {
                if is_fasta_aa_file(path) {
                    file_to_buffer(path)
                } else {
                    log::warn!(" encountering a not aa file : {:?}", path);
                    Vec::new()
                }
            }
        })
        .collect();
    //
    if log::log_enabled!(log::Level::Debug) {
        if cpu_start.is_ok() {
            let cpu_time = cpu_start.unwrap().try_elapsed();
            if cpu_time.is_ok() {
                log::debug!(
                    " files read! cpu_time (ms) : {:?}",
                    cpu_time.unwrap().as_millis()
                );
            }
        }
        let elapsed_t = start_t.elapsed().unwrap().as_millis() as f32;
        log::debug!(" files read! elapsed (ms) : {:?}", elapsed_t);
    }
    //
    // now we decompress and parse fasta buffers.
    //
    let mut to_be_processed: Vec<(&PathBuf, &[u8])> = Vec::with_capacity(pathb.len());
    let couples = pathb.iter().zip(&files_read);
    for c in couples {
        if !c.1.is_empty() {
            to_be_processed.push((c.0, c.1));
        }
    }
    // this is what we  wanted, fasta buffer parsing in // !!
    let seqseq = to_be_processed
        .par_iter()
        .map(|file| process_file(file.0, file.1).unwrap())
        .collect();
    // explicit drop to be able to monitor memory
    drop(to_be_processed);
    drop(files_read);
    //
    if log::log_enabled!(log::Level::Debug) {
        log::debug!("memory  : {:?}", memory_stats::memory_stats().unwrap());
        log::debug!(" end of process_files_group");
    }
    //
    seqseq
} // end of process_files_group

// recursive parallel processing of directores
pub(crate) fn process_dir_parallel_rec(
    state: &mut ProcessingState,
    datatype: &DataType,
    dir: &Path,
    filter_params: &FilterParams,
    block_size: usize,
    path_block: &mut Vec<PathBuf>,
    file_task: &(dyn Fn(&PathBuf, &[u8], &FilterParams) -> Vec<IdSeq> + Sync),
    sender: &crossbeam_channel::Sender<Vec<IdSeq>>,
) -> io::Result<usize> {
    //
    // We will process files in parallel by blocks of size block_size to control memory
    // and to balance cpu and memory with threads doing hnsw work
    // we want to keep track of how many file are processed.
    let mut nb_sent: usize = 0;
    let mut nb_entries = 0usize;
    //
    if std::fs::metadata(dir).is_err() {
        log::error!("\n\n process_dir_parallel_rec , dir : {:?}", dir);
        log::error!("\n directory do not exist or you do not have access permissions \n");
        std::panic!();
    }
    log::debug!("process_dir_parallel_rec , dir : {:?}", dir);
    //
    for entry in std::fs::read_dir(dir)? {
        if entry.is_err() {
            log::error!(
                "error reading entry in director : {:?}, entry num : {}",
                dir,
                nb_entries
            );
        }
        let entry = entry?;
        let pathb = entry.path();
        if pathb.is_dir() {
            nb_sent = process_dir_parallel_rec(
                state,
                datatype,
                &pathb,
                filter_params,
                block_size,
                path_block,
                file_task,
                sender,
            )?;
        } else {
            nb_entries += 1;
            // we have path of a file
            match pathb.metadata() {
                Ok(meta) => {
                    let f_len = meta.len();
                    log::trace!(
                        " filename : {:?}, length = {}",
                        pathb.file_name().unwrap_or_default(),
                        f_len
                    );
                }
                Err(_) => {
                    log::error!("cannot access to metadata of path : {:?}", pathb);
                }
            }
            if path_block.len() < block_size {
                // we push the file
                path_block.push(pathb.clone());
                // now if buffer is full we do the work, process send and empty buffer
                if path_block.len() == block_size {
                    let seqs = process_files_group(datatype, filter_params, path_block, file_task);
                    for mut seqfile in seqs {
                        for i in 0..seqfile.len() {
                            seqfile[i].rank = state.nb_seq;
                            seqfile[i].filerank = state.nb_file;
                            state.nb_seq += 1;
                        }
                        state.nb_file += 1;
                        nb_sent += 1;
                        sender.send(seqfile).unwrap();
                    }
                    if log::log_enabled!(log::Level::Info) && state.nb_file % 1000 == 0 {
                        log::info!(
                            "nb file processed : {}, nb sequences processed : {}",
                            state.nb_file,
                            state.nb_seq
                        );
                    }
                    if state.nb_file % 1000 == 0 {
                        println!(
                            "nb file processed : {}, nb sequences processed : {}",
                            state.nb_file, state.nb_seq
                        );
                    }
                    path_block.clear();
                }
            } else if path_block.len() == block_size {
                // cannot happen beccause we checked if we needed to send mesg and flush buffer
                std::panic!("process_dir_parallel buffer overflow , should not occur")
            }
        }
    } // end of for on entry
    // in fact in gtdb we have only one pure file in bottom directories!!! so this is useful only to realize it or in other cases.
    log::debug!(
        "process_dir_parallel got nb entries dir : {:?} nb_entries : {}",
        dir,
        nb_entries
    );
    //
    Ok(nb_sent)
} // end of process_dir_parallel_rec

// driver for parallel io processing.
// It calls the recursive function process_dir_parallel_rec and then treat residue of message not sent as buffer to send is not full
// at end of directory tree exploration
pub(crate) fn process_dir_parallel(
    state: &mut ProcessingState,
    datatype: &DataType,
    dirpath: &Path,
    filter_params: &FilterParams,
    nb_files_by_group: usize,
    file_task: &(dyn Fn(&PathBuf, &[u8], &FilterParams) -> Vec<IdSeq> + Sync),
    sender: &crossbeam_channel::Sender<Vec<IdSeq>>,
) -> io::Result<usize> {
    //
    log::info!(
        "files::process_dir_parallel : calling process_dir_parallel, nb_files in parallel : {}",
        nb_files_by_group
    );
    //
    let nb_files_at_entry = state.nb_file;
    let mut nb_sent_parallel_rec;
    let mut path_block = Vec::<PathBuf>::with_capacity(nb_files_by_group);
    //
    let res_nb_sent_parallel = process_dir_parallel_rec(
        state,
        datatype,
        dirpath,
        filter_params,
        nb_files_by_group,
        &mut path_block,
        &file_task,
        sender,
    );
    // we must treat residue in path_block if any
    match res_nb_sent_parallel {
        Ok(nb) => {
            nb_sent_parallel_rec = nb;
        }
        _ => {
            log::error!("\n some error occurred in process_dir_parallel_rec");
            std::panic!("\n some error occurred in process_dir_parallel_rec");
        }
    };
    // now must treat residue in path_block
    if !path_block.is_empty() {
        log::info!(
            "files::process_dir_parallel sending residue, size : {}",
            path_block.len()
        );
        let seqs = process_files_group(datatype, filter_params, &path_block, &file_task);
        for mut seqfile in seqs {
            for ids in &mut seqfile {
                ids.rank = state.nb_seq;
                ids.filerank = state.nb_file;
                state.nb_seq += 1;
            }
            state.nb_file += 1;
            nb_sent_parallel_rec += 1;
            sender.send(seqfile).unwrap();
            log::trace!(
                "sketchandstore_dir_compressedkmer parallel io case nb msg sent : {}",
                nb_sent_parallel_rec
            );
        }
        path_block.clear();
        if log::log_enabled!(log::Level::Info) && state.nb_file % 1000 == 0 {
            log::info!(
                "nb file processed : {}, nb sequences processed : {}",
                state.nb_file,
                state.nb_seq
            );
        }
        if state.nb_file % 1000 == 0 {
            println!(
                "nb file processed : {}, nb sequences processed : {}",
                state.nb_file, state.nb_seq
            );
        }
    }
    let nb_sent_parallel = state.nb_file - nb_files_at_entry;
    log::debug!(
        "process_dir_parallel parallel io case : all messages sent :{}",
        nb_sent_parallel_rec
    );

    //
    Ok(nb_sent_parallel)
} // end of process_dir_parallel
