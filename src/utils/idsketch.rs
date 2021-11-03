//! The module gathers structures used to identify sequences stroed in Hnsw (with probminhash hash data) from a id

use serde::{Deserialize, Serialize};
use serde_json::{to_writer};

use std::path::{PathBuf};
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter};
use std::time::{SystemTime};

use kmerutils::base::{sequence::*};


/// This structure is a summary of struct IdSeq
/// It is sent to sketcher thread, serialized and dump in a file so that its rank in dump file
/// is the data id used in Hnsw. So we can go back from a neighbour returned by hnsw to identity of sequence.
#[derive(Serialize,Deserialize, Clone)]
pub struct Id {
    /// path where seq was
    path : String,
    /// fasta id 
    fasta_id : String,
}

impl Id {
    pub fn new(path : &String, id : &String) -> Self {
        Id{ path : path.clone(), fasta_id : id.clone()}
    }

    /// get file path 
    pub fn get_path(&self) -> &String {
        &self.path
    }
    
    /// get sequence id in fasta file
    pub fn get_fasta_id(&self) -> &String {
        return &self.fasta_id;
    }
} // end of impl Id



/// 
/// This structure is used for returning info from function process_file
/// It stores all info on treated sequences 
/// 
pub struct IdSeq {
    /// as read is sequential we can identify uniquely sequence in hnsw
    pub(crate) rank : usize,
    /// But we do not know in which order files are read, so we store filename
    path : String,
    /// fasta id of genome Sketched as read in head of fasta record.
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

    /// get sequence length
    pub fn get_seq_len(&self) -> usize {
        self.seq.size()
    }
} // end of impl IdSea


/// We maintain an association of sequence Id and its length (useful to build a test?, anyway not a big cost )
#[derive(Clone, Serialize, Deserialize)]
pub struct ItemDict {
    /// sequence id
    id : Id, 
    /// sequence length
    len : usize,
}

impl ItemDict {
    ///
    pub fn new(id : Id, len : usize) -> Self {
        ItemDict{id, len}
    }
    ///
    pub fn get_id(&self) -> &Id { 
            &self.id
    }
    ///
    pub fn get_len(&self) -> usize {
         self.len
    }
}  // end of impl seqdict.0[n.d_id].get_id()


/// To keep track of sequence id by their rank.  
/// The rank of a sequence is its dataId in the Hnsw structure. 
/// So we can retrieve Seq description from Hnsw and the dictionary. 
pub struct SeqDict(pub Vec<ItemDict>);



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
    pub fn reload(filename : &String) -> Result<SeqDict, String>  {
        //
        log::info!("reloading database ids from : {}", filename);
        let start_t = SystemTime::now();
        //
        println!("reloading database ids from : {}", filename);
        let mut sequences =  Vec::<ItemDict>::with_capacity(100000);
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
        let stream = serde_json::Deserializer::from_reader(reader).into_iter::<ItemDict>();
        for value in stream {
            if value.is_err() {
                println!("nb id loaded {}", nbloaded);
                log::error!("nb id loaded {}", nbloaded);
                return Err(String::from("parsing error"));
            }
            sequences.push(value.unwrap());
            if log::log_enabled!(log::Level::Debug) && nbloaded <= 3 {
                log::debug!(" nbloaded : {:?}, sesqid path : ({}, {}), len : {}", nbloaded, sequences[nbloaded].get_id().path,
                        sequences[nbloaded].get_id().fasta_id, sequences[nbloaded].get_len());
            }
            nbloaded += 1;
        } 
        //
        log::info!("SeqDict, reloaded nb sequences : {:?}", nbloaded);
        let elapsed_t = start_t.elapsed().unwrap().as_secs() as f32;
        log::info!("SeqDict::reload : elapsed system time(s) {}", elapsed_t);
        //        
        return Ok(SeqDict{0:sequences});
    } // end of reload

    /// computes totol number of bases of database.
    pub fn get_total_length(&self) -> usize {
       self.0.iter().map(|s| s.len).sum()
    }

}  // end of impl SeqDict

