//! files utilities for dna processing
//!

use std::path::PathBuf;

use crate::utils::idsketch::*;
use crate::utils::parameters::*;
use kmerutils::base::sequence::*;

#[inline]
/// clones the sequence filtering out non ATCG
pub fn filter_out_n(seq: &[u8]) -> Vec<u8> {
    let mut filtered = Vec::<u8>::with_capacity(seq.len());

    for c in seq {
        let ch_up = (*c as char).to_ascii_uppercase();
        if ['A', 'C', 'T', 'G'].contains(&ch_up) {
            filtered.push(ch_up as u8);   // push the upper-case byte
        }
    }

    if log::log_enabled!(log::Level::Trace) && filtered.len() < seq.len() {
        let nb_n = seq.len() - filtered.len();
        log::trace!(
            "filtered nb non ACTG {}, fraction  {:1.3e}",
            nb_n,
            nb_n as f32 / seq.len() as f32
        );
    }
    filtered
} // end of filter_out_n

/* // TODO
 group process_file and process_dir in a structure that would maintain number of processed file
 This structure would maintain a triplet association (filename, rank in file, seqid)
 So it could be easy to have the sequence knowing the triplet
*/

/// opens parse fna files with needletail
/// extracts records , filters out capsid and send sequences to function process_dir to execute file_task to produce sequence
/// for any client.
/// This function returns as many IdSeq as there are sequences in file
pub fn process_file_by_sequence(pathb: &PathBuf, filter_params: &FilterParams) -> Vec<IdSeq> {
    let mut to_sketch = Vec::<IdSeq>::new();
    //
    log::debug!("processing file {}", pathb.to_str().unwrap());
    let mut nb_bases_file = 0;
    let mut nb_bases_encoded = 0;
    let mut nb_record: usize = 0;
    let alphabet2b = Alphabet2b::new();
    //
    let mut reader = needletail::parse_fastx_file(pathb).expect("expecting valid filename");
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", pathb.file_name().unwrap());
            std::process::exit(1);
        }
        nb_record += 1;
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        let file_seq = seqrec.seq();
        // we check for length, it happens (for aa file at least) that some sequence may be null
        if file_seq.len() > 0 {
            // process sequence if not capsid and not filtered out
            if !strid.contains("capsid") && !filter_params.filter(&file_seq) {
                let nb_bases = file_seq.len();
                nb_bases_file += nb_bases;
                let mut new_seq = Sequence::with_capacity(2, nb_bases);
                new_seq.encode_and_add(&file_seq, &alphabet2b);
                nb_bases_encoded += new_seq.size();
                // we have DNA seq for now
                // we encode path only fo the first seq. For files with hundreds millions of seq , memory
                let nullstr = String::from(""); // we will sketch the whole and loose id, so we spare memory
                let path = if to_sketch.is_empty() {
                    pathb.to_str().unwrap().to_string()
                } else {
                    nullstr.clone()
                };
                // we will sketch the whole and loose fasta id, so we spare memory
                let seqwithid = IdSeq::new(path, nullstr, SequenceType::SequenceDNA(new_seq));
                to_sketch.push(seqwithid);
            }
        } else {
            log::error!(
                "sequence of null length, file is {:?}, record num : {}",
                pathb,
                nb_record
            );
        }
        //
    } // end while
      //
    log::debug!(
        "process_file_by_sequence,  file : {:?}, nb_bases_encoded : {}, nb_record : {}",
        pathb,
        nb_bases_encoded,
        nb_record
    );
    log::debug!("nb total bases in file (with non ATCG) : {}", nb_bases_file);
    assert!(nb_bases_file >= nb_bases_encoded);
    // we must send to_sketch to some sketcher
    to_sketch
} // end of process_file_by_sequence

/// This function will parse with needletail (and do the decompressing).
/// We encode sequences in 2 bits the whole sequence on the fly; record after record; the whole bytes of file pathb contained in  bufread.
/// In this way process_buffer_by_sequence do not have any IO to do and can be called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq.  
/// The return of this function is a vector of size the number of record in the file.
pub fn process_buffer_by_sequence(
    pathb: &PathBuf,
    bufread: &[u8],
    filter_params: &FilterParams,
) -> Vec<IdSeq> {
    //
    log::debug!("process_buffer_by_sequence , file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    let alphabet2b = Alphabet2b::new();
    let mut nb_bases_file = 0;
    let mut nb_bases_encoded = 0;
    let mut nb_record: usize = 0;
    //
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    let mut record_num: u64 = 0;
    //
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!("got bd record in file {:?}", pathb.file_name().unwrap());
            std::process::exit(1);
        }
        nb_record += 1;
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        let file_seq = seqrec.seq();
        // we check for length it happens (for aa file at least) that some sequence may be null
        if file_seq.len() > 0 {
            // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
            if !strid.contains("capsid") && !filter_params.filter(&file_seq) {
                let nb_bases = file_seq.len();
                nb_bases_file += nb_bases;
                let mut new_seq = Sequence::with_capacity(2, nb_bases);
                new_seq.encode_and_add(&file_seq, &alphabet2b);
                // some checks , we filter non ACGT so length is less than record
                assert!(new_seq.size() <= nb_bases);
                // we have DNA seq for now
                nb_bases_encoded += new_seq.size();
                // we encode path only fo the first seq. For files with hundreds millions of seq , memory
                let nullstr = String::from("");
                let path = if to_sketch.is_empty() {
                    pathb.to_str().unwrap().to_string()
                } else {
                    nullstr.clone()
                };
                // we will sketch the whole and loose fasta id, so we spare memory
                let seqwithid = IdSeq::new(path, nullstr, SequenceType::SequenceDNA(new_seq));
                to_sketch.push(seqwithid);
            }
        } else {
            log::error!(
                "sequence of null length, file is {:?}, record num : {}",
                pathb,
                record_num
            );
        }
        record_num += 1;
    } // end while
      //
    log::debug!("process_buffer_by_sequence file : {:?}, nb_bases_file : {} nb_bases_encoded : {}, nb_record : {}", 
                        pathb, nb_bases_file, nb_bases_encoded, nb_record);
    log::debug!("nb total bases in file (with non ATCG) : {}", nb_bases_file);
    assert!(nb_bases_file >= nb_bases_encoded);
    //
    if log::log_enabled!(log::Level::Trace) {
        log::trace!("process_file, nb_sketched {} ", to_sketch.len());
    }
    //
    to_sketch
} // end of process_buffer_by_sequence

/// opens parse fna files with needletail extracts records , filters out capsid , encode in 2 bits the whome sequence on the fly,
/// record after record; the whole bytes of file pathb into one large sequence.
/// and send sequenceto function process_dir to execute file_task to produce sequence
/// for any client
/// The vector returned has size 1 as the sequence is concatenated
pub fn process_file_in_one_block(pathb: &PathBuf, filter_params: &FilterParams) -> Vec<IdSeq> {
    //
    log::debug!("process_file_in_one_block , file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    let mut nb_bases_file = 0;
    //
    let metadata = std::fs::metadata(pathb.clone());
    let nb_bases: usize;
    match metadata {
        Ok(metadata) => {
            if pathb.extension().is_some() && (pathb.extension().unwrap() == "gz") {
                // The decompressed file  is larger than the compressed one, expecting a factor 4 compression
                log::debug!("compressed file : {:?}", pathb);
                nb_bases = 4 * metadata.len() as usize;
            } else {
                log::debug!("uncompressed file : {:?}", pathb);
                nb_bases = metadata.len() as usize;
            };
        }
        Err(_) => {
            println!(
                "process_file_in_one_blprocess_dirock could not get length of file {} ",
                pathb.to_str().unwrap()
            );
            nb_bases = 1_000_000_000;
        }
    }
    log::trace!("processing file {}", pathb.to_str().unwrap());
    let bufread = std::io::BufReader::with_capacity(5_000_000, std::fs::File::open(pathb).unwrap());
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
        let file_seq = seqrec.seq();
        // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
        let _filter = filter_params.filter(&file_seq);
        if !strid.contains("capsid") {
            nb_bases_file += file_seq.len();
            new_seq.encode_and_add(&file_seq, &alphabet2b);
            // recall rank is set in process_dir beccause we should a have struct gatheing the 2 functions process_dir and process_file
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    log::debug!(
        "process_file_in_one_block file : {:?}, nb_bases : {}, seq len {}",
        pathb,
        nb_bases_file,
        new_seq.size()
    );
    log::debug!("nb total bases in file (with non ATCG) : {}", nb_bases_file);
    // should have equality except for non ACTG
    assert!(nb_bases_file >= new_seq.size());
    // we are at end of file, we have one large sequence for the whole file
    // we have DNA seq for now
    let seqwithid = IdSeq::new(
        pathb.to_str().unwrap().to_string(),
        String::from("-total-sequence"),
        SequenceType::SequenceDNA(new_seq),
    );
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    to_sketch
} // end of process_file_in_one_block

/// This function will parse with needletail (and do the decompressing).
/// We encode in 2 bits the whole fly file on the fly into a concatenated sequence, record after record; the whole bytes of file pathb contained in  bufread.
/// In this way process_buffer_in_one_block do not have any IO to do and can be called // for fasta parsing without disk constraints.
/// We nevertheless needs pathb to fill in IdSeq
/// The return of this function is a vector of size 1 as the sequence is concatenated
pub fn process_buffer_in_one_block(
    pathb: &PathBuf,
    bufread: &[u8],
    filter_params: &FilterParams,
) -> Vec<IdSeq> {
    //
    log::trace!("process_buffer_in_one_block , file : {:?}", pathb);
    //
    let mut to_sketch = Vec::<IdSeq>::new();
    let mut nb_bases_file = 0;
    let mut nb_record: usize = 0;
    //
    let mut reader = needletail::parse_fastx_reader(bufread).expect("expecting valid filename");
    // We allocate one large block tht will contain the whole filtered genome.
    let nb_bases = if pathb.extension().is_some() && (pathb.extension().unwrap() == "gz") {
        // The decompressed file  is larger than the compressed one
        bufread.len() * 4
    } else {
        bufread.len()
    };
    log::trace!(
        "allocating seq for {:?}, estimated nb bases : {}",
        pathb,
        nb_bases
    );
    // We allocate one large block tht will contain the whole filtered genome and we know the sequence size it will produce
    // if file is not compressed f_len is a majorant (as N and capsid are excluded), if file is compressed a compression of 2 can be expected
    let mut new_seq = Sequence::with_capacity(2, nb_bases);
    let alphabet2b = Alphabet2b::new();
    //
    while let Some(record) = reader.next() {
        if record.is_err() {
            println!(
                "process_buffer_in_one_block : got bad record in buffer : {:?}",
                pathb
            );
            std::process::exit(1);
        }
        nb_record += 1;
        // do we keep record ? we must get its id
        let seqrec = record.expect("invalid record");
        let id = seqrec.id();
        let strid = String::from_utf8(Vec::from(id)).unwrap();
        let file_seq = seqrec.seq();
        // process sequence if not capsid and not filtered out, in block mode we do not filter any at the moment
        let _filter = filter_params.filter(&file_seq);
        if !strid.contains("capsid") {
            // Our Kmers are 2bits encoded so we need to be able to encode sequence in 2 bits, so there is
            // this hack,  causing reallocation. seqrec.seq is Cow so drain does not seem an option.
            nb_bases_file += file_seq.len();
            new_seq.encode_and_add(&file_seq, &alphabet2b);
            if log::log_enabled!(log::Level::Trace) {
                log::trace!("process_file, nb_sketched {} ", to_sketch.len());
            }
        }
    }
    //
    new_seq.shrink_to_fit();
    // we are at end of file, we have one large sequence for the whole file
    log::debug!("decompressed seq for file : {:?}, nb bases encoded : {}, estimated nb bases : {}, nb_record : {}", pathb.file_name().unwrap_or_default(), new_seq.size(), nb_bases, nb_record);
    log::debug!("nb total bases in file (with non ATCG) : {}", nb_bases_file);
    assert!(nb_bases_file >= new_seq.size());
    // we have DNA seq for now
    let seqwithid = IdSeq::new(
        pathb.to_str().unwrap().to_string(),
        String::from("-total-sequence"),
        SequenceType::SequenceDNA(new_seq),
    );
    to_sketch.push(seqwithid);
    // we must send to_sketch to some sketcher
    to_sketch
} // end of process_buffer_in_one_block
