//! The module gathers structures used to identify sequences stroed in Hnsw (with probminhash hash data) from a id

use serde::{Serialize, Deserialize};


/// magic number identifying beginning of SeqDict dump
const MAGIC_SEQDICT : u32 = 0x56E289;



/// to keep track of sequence id by their rank. 
/// So we can retrieve Seq description from Hnsw and the dictionary 
/// TODO to be serialized and dumped
#[derive(Serialize, Deserialize)]
pub struct SeqDict(pub Vec<String>);



impl SeqDict {
    pub fn new(size : usize) -> Self {
        return SeqDict { 0 : Vec::with_capacity(size)};
    }

    /// serialize and dump
    /// At head of file we insert a Magic
    pub fn dump(filename : String) {

    } // end of dump

    /// reload from dump, necessary to run as a service
    pub fn reload(filename : String) {

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
