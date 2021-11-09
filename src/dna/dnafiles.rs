//! files utilities for dna processing
//! 


use std::fs::{self, DirEntry};

use crate::utils::{idsketch::*};
use crate::utils::{parameters::*};
use kmerutils::base::{sequence::*};


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
}  // end of filter_out_n



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
            // filtering cause reallocation. As we encode in 2 bits we cannot get something that does not fit. seqrec.seq is Cow so drain does not seem an option.
            let filtered = filter_out_n(&seqrec.seq());
            drop(seqrec);
            // we have DNA seq for now
            let new_seq = Sequence::new(&filtered,2);
            let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), strid,SequenceType::SequenceDNA(new_seq));
            to_sketch.push(seqwithid);
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_by_sequence




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
    match metadata {
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
        // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
        let _filter = filter_params.filter(&seqrec.seq());
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
    // we have DNA seq for now
    let new_seq = Sequence::new(&one_block_seq,2);
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceDNA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_in_one_block




/// opens parse fna files with needletail
/// extracts records , concat in a whole block and split in equal parts.
/// The split is done in 10 segments (hard coded at present time)
/// filters out capsid and send sequences to function process_dir to execute file_task to produce sequences
/// for any client
pub fn process_file_concat_split(file : &DirEntry, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let pathb = file.path();
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let metadata = fs::metadata(pathb.clone());
    let f_len : usize;
    match metadata {
        Ok(metadata) => { f_len = metadata.len() as usize;
                          log::debug!("file length : {}", f_len);
                        }
        Err(_)       => {
            println!("process_file_concat_split could not get length of file {} ", pathb.to_str().unwrap());
            f_len = 1_000_000_000;
        }
    }
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
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
        // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
        let _filter = filter_params.filter(&seqrec.seq());
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
    // we are at end of file, we have one large sequence and now we split it in 10 // TODO to parametrize
    let nb_split = 10;
    let block_size = one_block_seq.len() / nb_split;
    for i in 0..nb_split {
        let block_begin = i * block_size;
        let block_end = one_block_seq.len().min((i+1) * block_size) - 1;
        // generate an id
        let id = format!("total seq split {}", i);
        let new_seq = Sequence::new(&one_block_seq[block_begin..block_end],2);
        let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), id, SequenceType::SequenceDNA(new_seq));
        to_sketch.push(seqwithid);
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_concat_split

