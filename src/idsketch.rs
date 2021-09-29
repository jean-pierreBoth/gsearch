//! The module gathers structures used to identify sequences stroed in Hnsw (with probminhash hash data) from a id


use serde_json::{to_writer};

use std::path::{PathBuf};
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter};

/// magic number identifying beginning of SeqDict dump
const MAGIC_SEQDICT : u32 = 0x56E289;



/// to keep track of sequence id by their rank. 
/// So we can retrieve Seq description from Hnsw and the dictionary 
/// TODO to be serialized and dumped
pub struct SeqDict(pub Vec<String>);



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
            let id_str = serde_json::to_value(v).unwrap();
            let _ = to_writer(&mut writer, &MAGIC_SEQDICT).unwrap();
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
        let mut sequences =  Vec::<String>::with_capacity(100000);
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
        let stream = serde_json::Deserializer::from_reader(reader).into_iter::<serde_json::Value>();
        for value in stream {
            sequences.push(value.unwrap().to_string());
            if log::log_enabled!(log::Level::Debug) && nbloaded <= 3 {
                log::debug!(" nbloaded : {:?}, sesqid : {:?}", nbloaded, sequences[nbloaded]);
            }
            nbloaded += 1;
        } 
        //
        log::info!("SeqDict, reloaded nb sequences : {:?}", nbloaded);
        //        
        return Ok(SeqDict{0:sequences});
    } // end of reload

}  // end of impl SeqDict


/// The structure to be embedded in Hnsw.
/// The distance declared in Hnsw will be be on IdSketch, operating on the field sig
/// The id is provided by reading thread and is stored separetely in a Dictionary to be serialized and dumped as the Hnsw structure.
pub struct IdSketch {
    /// (unique) id of genome Sketched
    id : usize,
    /// probminhash signature
    sig : Vec<u32>,
} // end of IdSketch
