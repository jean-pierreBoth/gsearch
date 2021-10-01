//! This file contains directory exploration and fasta file selection
#![allow(unused)]

use std::io;
use std::fs::{self, DirEntry};
use std::path::Path;
use kmerutils::base::{sequence::*};

use crate::idsketch::{IdSeq};


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
pub fn process_file(file : &DirEntry)  -> Vec<IdSeq> {
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
        if strid.find("capsid").is_none() {
            // Our Kmers are 2bits encoded so we need to be able to encode sequence in 2 bits, so there is 
            // this hack,  causing reallocation. seqrec.seq is Cow so drain does not seem an option.
            let filtered = filter_out_n(&seqrec.seq());
            let newseq = Sequence::new(&filtered, 2);
            if log::log_enabled!(log::Level::Trace) && filtered.len() < seqrec.seq().len() {
                let nb_n = seqrec.seq().len() - filtered.len();
                log::trace!("filtered nb non ACTG {}, fraction  {:1.3e}", seqrec.seq().len() - filtered.len(), nb_n as f32/seqrec.seq().len() as f32);
            }
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            let seqwithid = IdSeq::new(pathb.to_str().unwrap().to_string(), strid, newseq);
            to_sketch.push(seqwithid);
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    // we must send to_sketch to some sketcher
    return to_sketch;
} // end of process_file



/// scan directory recursively, executing function file_task on each file.
/// adapted from from crate fd_find
pub fn process_dir(dir: &Path, file_task: &dyn Fn(&DirEntry) -> Vec<IdSeq>, sender : &crossbeam_channel::Sender::<Vec<IdSeq>>) -> io::Result<usize> {
    let mut nb_seq_processed = 0;
    //
    // we checked that we have a directory
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            nb_seq_processed += process_dir(&path, file_task, sender)?;
        } else {
            // check if entry is a fasta.gz file
            if is_fasta_file(&entry) {
                let mut to_sketch = file_task(&entry);
                // put a rank id in sequences, now we have full information of where do the sequence come from
                for i in 0..to_sketch.len() {
                    to_sketch[i].rank = nb_seq_processed + i;
                }
                nb_seq_processed += to_sketch.len();
                // we must send to_sketch into channel to upper thread
                sender.send(to_sketch).unwrap();
            }
        }
    }
    if nb_seq_processed > 0 && nb_seq_processed % 50000 == 0 {
        log::info!("processed nb sequences : {}", nb_seq_processed);
    }
    log::trace!("processed nb sequences : {}", nb_seq_processed);
    //
    drop(sender);
    //
    Ok(nb_seq_processed)
}  // end of visit_dirs

