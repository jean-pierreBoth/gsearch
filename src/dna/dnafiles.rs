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
/// for any client.
/// This function returns as many IdSeq as there are sequences in file
pub fn process_file_by_sequence(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let alphabet2b = Alphabet2b::new();
    //
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
            let nb_bases = seqrec.seq().len();
            let mut new_seq = Sequence::with_capacity(2, nb_bases);
            new_seq.encode_and_add(&seqrec.seq(), &alphabet2b);
            drop(seqrec);
            // we have DNA seq for now
            let nullid = String::from(""); // we will sketch the whole and loose id, so we spare memory
            let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), nullid,SequenceType::SequenceDNA(new_seq));
            to_sketch.push(seqwithid);
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_by_sequence




/// This function will parse with needletail (and do the decompressing).
/// We encode sequences in 2 bits the whole sequence on the fly; record after record; the whole bytes of file pathb contained in  bufread.
/// In this way process_buffer_in_one_block do not have any IO to do and can be called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq.  
/// The return of this function is a vector of size the number of record in the file.
pub fn process_buffer_by_sequence(pathb : &PathBuf, bufread : &[u8], filter_params : &FilterParams)  -> Vec<IdSeq> {
    //
    log::debug!("process_buffer_by_sequence , file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    let alphabet2b = Alphabet2b::new();
    //
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
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
        if strid.find("capsid").is_none()  && !filter_params.filter(&seqrec.seq()) {
            let nb_bases = seqrec.seq().len();
            let mut new_seq = Sequence::with_capacity(2, nb_bases);
            new_seq.encode_and_add(&seqrec.seq(), &alphabet2b);
            drop(seqrec);
            // we have DNA seq for now
            let nullid = String::from(""); // we will sketch the whole and loose id, so we spare memory
            let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), nullid,SequenceType::SequenceDNA(new_seq));
            to_sketch.push(seqwithid);
        }
    }
    if log::log_enabled!(log::Level::Trace) {
        log::trace!("process_file, nb_sketched {} ", to_sketch.len());
    }
    //
    return to_sketch;
} // end of process_buffer_by_sequence



/// opens parse fna files with needletail extracts records , filters out capsid , encode in 2 bits the whome sequence on the fly,
/// record after record; the whole bytes of file pathb into one large sequence.
/// and send sequenceto function process_dir to execute file_task to produce sequence
/// for any client
/// The vector returned has size 1 as the sequence is concatenated
pub fn process_file_in_one_block(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
    //
    log::debug!("process_file_in_one_block , file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let metadata = std::fs::metadata(pathb.clone());
    let nb_bases : usize;
    match metadata {
        Ok(metadata) => { 
            if pathb.extension().is_some() && (pathb.extension().unwrap() == "gz") {
                // The decompressed file  is larger than the compressed one, expecting a factor 4 compression
                log::debug!("compressed file : {:?}", pathb);
                nb_bases = 4 * metadata.len() as usize;
            }
            else {
                log::debug!("uncompressed file : {:?}", pathb);
                nb_bases = metadata.len() as usize;
            };
        },
        Err(_)       => {
            println!("process_file_in_one_blprocess_dirock could not get length of file {} ", pathb.to_str().unwrap());
            nb_bases = 1_000_000_000;
        }
    }
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let bufread = std::io::BufReader::with_capacity(5_000_000,std::fs::File::open(pathb).unwrap());
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome and we know the sequence size it will produce
    // if file is not compressed f_len is a majorant (as N and capsid are excluded), if file is compressed a compression of 2 can be expected
    let mut new_seq = Sequence::with_capacity(2, nb_bases);
    let alphabet2b = Alphabet2b::new();
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
            new_seq.encode_and_add(&seqrec.seq(), &alphabet2b);
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we are at end of file, we have one large sequence for the whole file
    // we have DNA seq for now
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceDNA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file_in_one_block





/// This function will parse with needletail (and do the decompressing).
/// We encode in 2 bits the whole fly file on the fly into a concatenated sequence, record after record; the whole bytes of file pathb contained in  bufread.
/// In this way process_buffer_in_one_block do not have any IO to do and can be called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq
/// The return of this function is a vector of size 1 as the sequence is concatenated
pub fn process_buffer_in_one_block(pathb : &PathBuf, bufread : &[u8], filter_params : &FilterParams)  -> Vec<IdSeq> {
    //
    log::trace!("process_buffer_in_one_block , file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
    let nb_bases;
    if pathb.extension().is_some() && (pathb.extension().unwrap() == "gz") {
        // The decompressed file  is larger than the compressed one
        nb_bases = bufread.len() * 4;
    }
    else {
        nb_bases = bufread.len();
    }
    log::trace!("allocating seq for {:?}, estimated nb bases : {}", pathb, nb_bases);
    // We allocate one large block tht will contain the whole filtered genome and we know the sequence size it will produce
    // if file is not compressed f_len is a majorant (as N and capsid are excluded), if file is compressed a compression of 2 can be expected
    let mut new_seq = Sequence::with_capacity(2, nb_bases);
    let alphabet2b = Alphabet2b::new();
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
            new_seq.encode_and_add(&seqrec.seq(), &alphabet2b);
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    new_seq.shrink_to_fit();
    // we are at end of file, we have one large sequence for the whole file
    log::debug!("decompressed seq for file : {:?}, nb bases : {}, estimated nb bases : {}", pathb.file_name().unwrap_or_default(), new_seq.size(), nb_bases);
    // we have DNA seq for now
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceDNA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    to_sketch
}  // end of process_buffer_in_one_block



