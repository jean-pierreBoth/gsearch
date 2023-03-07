//! constains answer to requests
//! 
//! 
use std::fs::{File};
use std::io::{Write,BufWriter};

use crate::utils::idsketch::{ItemDict, SeqDict};

use hnsw_rs::prelude::*;


/// An answer has a rank (in which it is processed), the fasta Id corresponding to request sequence, and
/// the asked list of neighbours.
/// The neighbours are identified by an id in the database. To retrieve the fasta identity 
/// of the neighbour we will need the SeqDict
pub struct ReqAnswer<'a> {
    rank : usize,
    /// request id
    req_item : ItemDict,
    ///
    neighbours : &'a Vec<Neighbour>,
}


impl <'a> ReqAnswer<'a> {
    pub fn new(rank : usize, req_item : ItemDict, neighbours : &'a Vec<Neighbour>) -> Self {
        ReqAnswer { rank, req_item, neighbours}
    }

    /// dump answers to a File. 
    /// We dump only answers with distance less than threshold to help visual synthesis of reult.
    /// Typically keep only distance less than 0.98 with kmer size=12 is sufficient to get rid of garbage.
    pub(crate) fn dump(&self, seqdict : &SeqDict, threshold : f32, out : &mut BufWriter<File>) -> std::io::Result<usize> {
        // dump rank , fasta_id
        let has_match = self.neighbours.iter().any(|&n| n.distance <= threshold);
        let mut nb_match = 0;
        if has_match {
            write!(out, "\n\n {} path {}, fasta_id {}, len : {}", self.rank, self.req_item.get_id().get_path(), 
                    self.req_item.get_id().get_fasta_id(), self.req_item.get_len())?;
            for n in self.neighbours {
                // get database identification of neighbour
                if n.distance  < threshold {
                    nb_match += 1;
                    let database_id = seqdict.0[n.d_id].get_id().get_path();
                    write!(out, "\n\t distance : {:.3E}  answer fasta id {}", n.distance, database_id)?;
                    log::debug!(" \t data id : {}", n.d_id);
                    write!(out, "\n\t\t answer fasta id {}, seq len : {}", seqdict.0[n.d_id].get_id().get_fasta_id(), seqdict.0[n.d_id].get_len())?;
                }
            }
        } // end match
        Ok(nb_match)
    } // end of dump


    pub fn get_request_id(&self) -> &ItemDict {
        &self.req_item
    } // end of get_item_dict


}   // end of ReqAnswer


