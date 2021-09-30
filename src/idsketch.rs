//! The module gathers structures used to identify sequences stroed in Hnsw (with probminhash hash data) from a id

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};

use std::path::{PathBuf};
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter};



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
        let filepath = PathBuf::from(filename);
        let fileres = OpenOptions::new().write(true).create(true).truncate(true).open(&filepath);
        if fileres.is_err() {
            log::error!("SeqDict dump : dump could not open file {:?}", filepath.as_os_str());
            println!("SeqDict dump: could not open file {:?}", filepath.as_os_str());
            return Err("SeqDDeserializerut writer bad magic").unwrap();
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

