//! This file contains directory exploration and fasta file selection
#![allow(unused)]

use std::io;
use std::fs::{self, DirEntry};
use std::path::Path;
use kmerutils::base::{sequence::*};

use super::idsketch::{IdSeq};

/// To keep track of processed file and sequence processed
pub struct ProcessingState {
    /// nb sequences processed
    nb_seq : usize,
    /// nb file processed
    nb_file : usize,
} 


impl ProcessingState {
    pub fn new() -> Self {
        ProcessingState{nb_seq : 0, nb_file : 0}
    } // end of new
     

} // end of ProcessingState





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




// returns true if file is a fasta file (possibly gzipped)
// filename are of type GCA[GCF]_000091165.1_genomic.fna.gz
pub fn is_fasta_file(file : &DirEntry) -> bool {
    let filename = file.file_name().into_string().unwrap();
    if filename.ends_with("fna.gz")|| filename.ends_with("fa.gz") || filename.ends_with("fasta.gz") {
        return true;
    }
    else { 
        return false;
    }
}  // end of is_fasta_file


#[inline]
/// clones the sequence filtering out non ATCG
pub fn filter_out_n(seq : &[u8]) -> Vec<u8> {
    let mut filtered = Vec::<u8>::with_capacity(seq.len());
    for c in seq {
        if ['A','C', 'T', 'G'].contains(&(*c as char))  {
            filtered.push(*c);
        }
    }
    if log::log_enabled!(log::Level::Trace) && filtered.len() < seq.len() {
        let nb_n = seq.len() - filtered.len();
        log::trace!("filtered nb non ACTG {}, fraction  {:1.3e}", nb_n , nb_n as f32/seq.len() as f32);
    }
    return filtered;
}  // end of keep_atcg



/* // TODO
 group process_file and process_dir in a structure that would maintain number of processed file
 This structure would maintain a triplet association (filename, rank in file, seqid) 
 So it could be easy to have the sequence knowing the triplet
*/


/// opens parse fna files with needletail
/// extracts records , filters out capsid and send sequences to function process_dir to execute file_task to produce sequence
/// for any client
pub fn process_file_by_sequence(file : &DirEntry, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let pathb = file.path();
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", file.file_name());
            std::process::exit(1);
        }
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        // process sequence if not capsid and not filtered out
        if strid.find("capsid").is_none() && !filter_params.filter(&seqrec.seq()) {
            // Our Kmers are 2bits encoded so we need to be able to encode sequence in 2 bits, so there is 
            // this hack,  causing reallocation. seqrec.seq is Cow so drain does not seem an option.
            let old_len = seqrec.seq().len();
            let filtered = filter_out_n(&seqrec.seq());
            drop(seqrec);
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), strid, Sequence::new(&filtered,2));
            to_sketch.push(seqwithid);
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file


/// opens parse fna files with needletail
/// extracts records , filters out capsid and send sequences to function process_dir to execute file_task to produce sequence
/// for any client
pub fn process_file_in_one_block(file : &DirEntry, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let pathb = file.path();
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let metadata = fs::metadata(pathb.clone());
    let f_len : usize;
    match(metadata) {
        Ok(metadata) => { f_len = metadata.len() as usize;
                          log::debug!("file length : {}", f_len);
                        }
        Err(_)       => {
            println!("process_file_in_one_blprocess_dirock could not get length of file {} ", pathb.to_str().unwrap());
            f_len = 1_000_000_000;
        }
    }
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
    // TODO We should get the file length to optimize the length
    let mut one_block_seq = Vec::<u8>::with_capacity(f_len as usize);
    //
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", file.file_name());
            std::process::exit(1);
        }
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        // process sequence if not capsid and not filtered out
        if strid.find("capsid").is_none() {
            // Our Kmers are 2bits encoded so we need to be able to encode sequence in 2 bits, so there is 
            // this hack,  causing reallocation. seqrec.seq is Cow so drain does not seem an option.
            one_block_seq.append(&mut filter_out_n(&seqrec.seq()));
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we are at end of file, we have one large sequence for the whole file
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), Sequence::new(&one_block_seq,2));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file



/// scan directory recursively, executing function file_task on each file.
/// adapted from from crate fd_find
pub fn process_dir(state : &mut ProcessingState, dir: &Path, filter_params : &FilterParams, file_task: &dyn Fn(&DirEntry, &FilterParams) -> Vec<IdSeq>, sender : &crossbeam_channel::Sender::<Vec<IdSeq>>) -> io::Result<usize> {
    let mut nb_seq_processed = 0;
    let mut nb_file_processed = 0;
    //
    // we checked that we have a directory
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            nb_seq_processed += process_dir(state, &path, filter_params, file_task, sender)?;
        } else {
            // check if entry is a fasta.gz file
            if is_fasta_file(&entry) {
                let mut to_sketch = file_task(&entry, filter_params);
                // put a rank id in sequences, now we have full information of where do the sequence come from
                for i in 0..to_sketch.len() {
                    state.nb_seq += 1;
                    to_sketch[i].rank = state.nb_seq;
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
    Ok(nb_seq_processed)
}  // end of visit_dirs

