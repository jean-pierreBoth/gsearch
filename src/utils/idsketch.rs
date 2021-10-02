//! The module gathers structures used to identify sequences stroed in Hnsw (with probminhash hash data) from a id
#![allow(unused)]

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};

use std::path::{PathBuf};
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter};

use kmerutils::base::{sequence::*};


/// This structure is a summary of struct IdSeq
/// It is sent to sketcher thread, serialized and dump in a file so that its rank in dump file
/// is the data id used in Hnsw. So we can go back from a neighbour returned by hnsw to identity of sequence.
#[derive(Serialize,Deserialize)]
pub struct Id {
    /// path where seq was
    path : String,
    ///
    fasta_id : String,
}

impl Id {
    pub fn new(path : &String, id : &String) -> Self {
        Id{ path : path.clone(), fasta_id : id.clone()}
    }

} // end of impl Id



/// parameters used in constructing database and for requests.
/// The same values must be used for construction and requests so the structure is json serialized
#[derive(Serialize,Deserialize)]
pub struct SketcherParams {
    kmer_size : usize,
    sketch_size : usize,
}

impl SketcherParams {
    /// allocator
    pub fn new(kmer_size : usize, sketch_size : usize) -> Self {
        SketcherParams{kmer_size,sketch_size}
    }
    /// returns kmer size
    pub fn get_kmer_size(&self) -> usize {
        self.kmer_size
    }

    /// return sketch size
    pub fn get_sketch_size(&self) -> usize {
        self.sketch_size
    }  
    
    /// serialized dump
    pub fn dump_json(&self, filename : &String) -> Result<(), String> {
        //
        let filepath = PathBuf::from(filename.clone());
        //
        log::info!("dumping sketching parameters in json file : {}", filename);
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SketcherParams dump : dump could not open file {:?}", filepath.as_os_str());
            println!("SketcherParams dump: could not open file {:?}", filepath.as_os_str());
            return Err("SketcherParams dump failed".to_string());
        }
        // 
        let mut writer = BufWriter::new(fileres.unwrap());
        let _ = to_writer(&mut writer, &self).unwrap();
        //
        Ok(())
    } // end of dump


    /// reload from a json dump
    pub fn reload_json(filename : String) -> Result<SketcherParams, String> {
        let filepath = PathBuf::from(filename);
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SketcherParams reload_json : reload could not open file {:?}", filepath.as_os_str());
            println!("SketcherParams reload_json: could not open file {:?}", filepath.as_os_str());
            return Err("SketcherParams reload_json could not open file".to_string());            
        }
        //
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        let sketch_params:SketcherParams = serde_json::from_reader(reader).unwrap();      
        //
        Err("not yet implemented".to_string())
    } // end of reload_json


} // end of impl SketcherParams



/// 
/// This structure is used for returning info from function process_file
/// It stores all info on sequence. 
/// 
pub struct IdSeq {
    /// as read is sequential we can identify uniquely sequence in hnsw
    pub(crate) rank : usize,
    /// But we do not know in which order files are read, so we strore filename
    path : String,
    /// id of genome Sketched as read in head of fasta record.
    id : String,
    /// Sequence compressed to 2 bit / base
    seq : Sequence
}  // end of IdSeq


impl IdSeq {
    ///
    pub fn new(path : String, id: String, seq : Sequence) -> Self {
        IdSeq{rank : 0, path, id, seq}
    }
    /// get file path 
    pub fn get_path(&self) -> &String {
        &self.path
    }
    
    /// get fasta id
    pub fn get_fasta_id(&self) -> &String {
        &self.id
    }

    pub fn get_rank(&self) -> usize {
        self.rank
    }

    pub fn get_sequence(&self) -> &Sequence {
        &self.seq
    }
} // end of impl IdSea


/// to keep track of sequence id by their rank. 
/// So we can retrieve Seq description from Hnsw and the dictionary 
/// 
pub struct SeqDict(pub Vec<Id>);



impl SeqDict {
    pub fn new(size : usize) -> Self {
        return SeqDict { 0 : Vec::with_capacity(size)};
    }

    /// serialize and dump
    /// At head of file we insert a Magic MAGIC_SEQDICT
    pub fn dump(&self, filename : String) -> Result<(), String> {
        let filepath = PathBuf::from(filename.clone());
        //
        log::info!("dumping sequence dicitionary in json file : {}", filename);
        //
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SeqDict dump : dump could not open file {:?}", filepath.as_os_str());
            println!("SeqDict dump: could not open file {:?}", filepath.as_os_str());
            return Err("SeqDict Deserializer dump failed").unwrap();
        }
        let mut writer = BufWriter::new(fileres.unwrap());
        for v in &self.0 {
            // v is and Id struct
            let _ = to_writer(&mut writer, &v).unwrap();
        }
        //
        println!("SeqDict dumped sequence dictionary, nb seq : {}", self.0.len());
        //
        return Ok(());
    } // end of dump


    /// reload from dump to avoid parsing again files. 
    /// To be used with reload of Hnsw structure to run as a service
    pub fn reload(filename : String) -> Result<SeqDict, String>  {
        //
        let mut sequences =  Vec::<Id>::with_capacity(100000);
        let filepath = PathBuf::from(filename);
        let fileres = OpenOptions::new().read(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SeqDict reload : reload could not open file {:?}", filepath.as_os_str());
            println!("SeqDict reload: could not open file {:?}", filepath.as_os_str());
            return Err("SeqDict reload could not open file".to_string());            
        }
        //
        let mut nbloaded = 0;
        let loadfile = fileres.unwrap();
        let reader = BufReader::new(loadfile);
        // we must deserialize Id structs 
        let stream = serde_json::Deserializer::from_reader(reader).into_iter::<Id>();
        for value in stream {
            sequences.push(value.unwrap());
            if log::log_enabled!(log::Level::Debug) && nbloaded <= 3 {
                log::debug!(" nbloaded : {:?}, sesqid path : ({}, {})", nbloaded, sequences[nbloaded].path,sequences[nbloaded].fasta_id);
            }
            nbloaded += 1;
        } 
        //
        log::info!("SeqDict, reloaded nb sequences : {:?}", nbloaded);
        //        
        return Ok(SeqDict{0:sequences});
    } // end of reload

}  // end of impl SeqDict

