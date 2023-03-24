//! files utilities for dna processing
//! 


use std::path::PathBuf;

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
pub fn process_file_by_sequence(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", pathb.file_name().unwrap());
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
pub fn process_file_in_one_block(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let metadata = std::fs::metadata(pathb.clone());
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
    let bufread = std::io::BufReader::with_capacity(5_000_000,std::fs::File::open(pathb).unwrap());
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
    // TODO We should get the file length to optimize the length
    let mut one_block_seq = Vec::<u8>::with_capacity(f_len as usize);
    //
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", pathb.file_name().unwrap());
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





/// This function will parse with needletail (and do the decompressing) the whole bytes of file pathb contained in  bufread.
/// In this way process_buffer_in_one_block do not have any IO to do and can b called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq
pub fn process_buffer_in_one_block(pathb : &PathBuf, bufread : &[u8], filter_params : &FilterParams)  -> Vec<IdSeq> {
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
    let mut one_block_seq : Vec::<u8>;
    if pathb.ends_with(".gz") {
        // The decompressed file  is larger than the compressed one
        one_block_seq = Vec::<u8>::with_capacity(bufread.len() * 2);
    }
    else {
        one_block_seq = Vec::<u8>::with_capacity(bufread.len());
    }
    //
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("process_buffer_in_one_block : got bad record in buffer : {:?}", pathb);
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
    log::debug!("decompressed seq for file : {:?}, len is : {}", pathb.file_name().unwrap_or_default(), one_block_seq.len());
    // we have DNA seq for now
    let new_seq = Sequence::new(&one_block_seq,2);
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceDNA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    to_sketch
}  // end of process_buffer_in_one_block



/// opens parse fna files with needletail
/// extracts records , concat in a whole block and split in equal parts.
/// The split is done in 10 segments (hard coded at present time)
/// filters out capsid and send sequences to function process_dir to execute file_task to produce sequences
/// for any client
pub fn process_file_concat_split(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let metadata = std::fs::metadata(pathb.clone());
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
            println!("got bd record in file {:?}", pathb.as_path());
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
    // we are at end of file, we have one large sequence and now we split it in 10M
    let seq_len =  one_block_seq.len();
    let max_size = 10_000_000;
    let nb_split = if seq_len % max_size == 0 { seq_len / max_size} else { 1 + seq_len / max_size};
    log::debug!("split seq length : {:#} in nb blocks : {:#}", seq_len, nb_split);
    let block_size = seq_len / nb_split;
    for i in 0..nb_split {
        let block_begin = i * block_size;
        let block_end = seq_len.min((i+1) * block_size) - 1;
        // generate an id
        let id = format!("total seq split {}", i);
        let new_seq = Sequence::new(&one_block_seq[block_begin..block_end],2);
        let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), id, SequenceType::SequenceDNA(new_seq));
        to_sketch.push(seqwithid);
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_concat_split

