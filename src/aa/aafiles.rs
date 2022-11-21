//! stuff related to rna files work


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
pub fn process_aafile_in_one_block(pathb : &PathBuf, filter_params : &FilterParams)  -> Vec<IdSeq> {
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
    let mut reader = needletail::parse_fastx_file(&pathb).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome. 
    // TODO We should get the file length to optimize the length
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
            one_block_seq.append(&mut filter_out_non_aa(&seqrec.seq()));
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we are at end of file, we have one large sequence for the whole file
    // we have DNA seq for now
    let new_seq = kmeraa::SequenceAA::new(&one_block_seq);
    let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), String::from("total sequence"), SequenceType::SequenceAA(new_seq));
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_aafile_in_one_block


