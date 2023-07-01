//! stuff related to amino acid (AA) files work


use std::path::PathBuf;

use kmerutils::aautils::kmeraa;

use crate::utils::{idsketch::*};
use crate::utils::{parameters::*};


/// filters out poential non AA letters. In fact it is the * that is really searched.
pub fn filter_out_non_aa(seq : &[u8]) -> Vec<u8> {
    let mut filtered = Vec::<u8>::with_capacity(seq.len());
    let alphabet = kmeraa::Alphabet::new();
    for c in seq {
        if alphabet.is_valid_base(*c) {
            filtered.push(*c);
        }
    }
    if log::log_enabled!(log::Level::Trace) && filtered.len() < seq.len() {
        let nb_n = seq.len() - filtered.len();
        log::trace!("filtered nb non AA letters {}, fraction  {:1.3e}", nb_n , nb_n as f32/seq.len() as f32);
    }
    return filtered;
}  // end of filter_out_non_aa



/// opens parse fna files with needletail
/// extracts records , filters out capsid and send sequences to function process_dir to execute file_task to produce sequence
/// for any client
pub(crate) fn process_aafile_in_one_block(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
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
            println!("process_aafile_in_one_block could not get length of file {} ", pathb.to_str().unwrap());
            f_len = 100_000_000;
        }
    }
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
    // TODO We should get the file length to optimize the length
    let mut one_block_seq : Vec::<u8>;
    if pathb.ends_with(".gz") {
        // The decompressed file  is larger than the compressed one
        one_block_seq = Vec::<u8>::with_capacity(f_len * 2);
    }
    else {
        one_block_seq = Vec::<u8>::with_capacity(f_len);
    }    
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
            one_block_seq.append(&mut filter_out_non_aa(&seqrec.seq()));
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we are at end of file, we have one large sequence for the whole file
    // we have AA seq for now
    let new_seq = kmeraa::SequenceAA::new(&one_block_seq);
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceAA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_aafile_in_one_block




/// This function will parse with needletail (and do the decompressing) the whole bytes of file pathb contained in  bufread.
/// In this way process_aabuffer_in_one_block do not have any IO to do and can b called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq
pub fn process_aabuffer_in_one_block(pathb : &PathBuf, bufread : &[u8], filter_params : &FilterParams)  -> Vec<IdSeq> {
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
            println!("process_buffer_in_one_block : got bd record in buffer");
            std::process::exit(1);
        }
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
        let _filter = filter_params.filter(&seqrec.seq());
        if strid.find("capsid").is_none() {
            one_block_seq.append(&mut filter_out_non_aa(&seqrec.seq()));
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we are at end of file, we have one large sequence for the whole file
    log::debug!("decompressed seq for file : {:?}, seq len is : {}, file length : {}", pathb.file_name().unwrap_or_default(), one_block_seq.len(), bufread.len());
    // we have DNA seq for now
    let new_seq = kmeraa::SequenceAA::new(&one_block_seq);
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceAA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    to_sketch
}  // end of process_buffer_in_one_block



/// This function will parse with needletail (and do the decompressing).
/// We encode sequences in 2 bits the whole sequence on the fly; record after record; the whole bytes of file pathb contained in  bufread.
/// In this way process_buffer_by_sequence do not have any IO to do and can be called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq.  
/// The return of this function is a vector of size the number of record in the file.
pub fn process_aabuffer_by_sequence(pathb : &PathBuf, bufread : &[u8], filter_params : &FilterParams)  -> Vec<IdSeq> {
    //
    log::debug!("process_buffer_by_sequence , aa file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    let alphabet = kmeraa::Alphabet::new();
    //
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    let mut record_num: u64 = 0;
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
        // it happens that aa files have null sequences!!
        if seqrec.seq().len() > 0 {
            if strid.find("capsid").is_none()  && !filter_params.filter(&seqrec.seq()) {
                // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
                let new_seq = kmeraa::SequenceAA::new_filtered(&seqrec.seq(), &alphabet);
                if new_seq.len() > 0 {
                    let nullid = String::from(""); // we will sketch the whole and loose id, so we spare memory
                    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), nullid,SequenceType::SequenceAA(new_seq));
                    to_sketch.push(seqwithid);
                }
                else {
                    log::warn!("null encoded sequence in file : {:?}, record num : {}, record seq : {:?}", pathb, record_num, std::str::from_utf8(&seqrec.seq()).unwrap());
                    log::warn!("sequence id : {}", strid);
                }
            }
        }
        else {
            log::error!("sequence of null length, file is {:?}, record num : {}", pathb, record_num);
        }
        //
        drop(seqrec);
        record_num += 1;
    }
    if log::log_enabled!(log::Level::Trace) {
        log::trace!("process_file, nb_sketched {} ", to_sketch.len());
    }
    //
    return to_sketch;
} // end of process_buffer_by_sequence


/// opens parse fna files with needletail
/// extracts records , filters out capsid and send sequences to function process_dir to execute file_task to produce sequence
/// for any client.
/// This function returns as many IdSeq as there are sequences in file
pub fn process_aafile_by_sequence(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    log::debug!("processing file {}", pathb.to_str().unwrap());
    let alphabet = kmeraa::Alphabet::new();
    //
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    let mut record_num: u64 = 0;
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
        // it happens that aa files have null sequences!!
        if seqrec.seq().len() > 0 {
            // process sequence if not capsid and not filtered out
            if strid.find("capsid").is_none() && !filter_params.filter(&seqrec.seq()) {
                let new_seq = kmeraa::SequenceAA::new_filtered(&seqrec.seq(), &alphabet);
                // we have AA seq for now
                if new_seq.len() > 0 {
                    let nullid = String::from(""); // we will sketch the whole and loose id, so we spare memory
                    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), nullid,SequenceType::SequenceAA(new_seq));
                    to_sketch.push(seqwithid);
                }
                else {
                    log::warn!("null encoded sequence in file : {:?}, record num : {}, record seq : {:?}", pathb, record_num, std::str::from_utf8(&seqrec.seq()).unwrap());
                    log::warn!("sequence id : {}", strid);
                }
            }
        }
        else {
            log::error!("sequence of null length, file is {:?}, record num : {}", pathb, record_num);
        }
        //
        drop(seqrec);
        record_num += 1;
    } // end while 
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_by_sequence

