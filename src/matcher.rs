//! A module to quantify match between sequence.  
//! 
//! For each file/genome we maintain a list of matches with other files/genomes
//! A match between sequences is a couple of sequence one from request genome, one from the database
//! for which we found a distance below some significant threshold.
//! 
//! A match between genomes is a list of matches between couple of sequences, one of the request genome, 
//! the other from the database genomes.

#![allow(unused)]

use std::collections::HashMap;

use crate::utils::*;

#[derive(Clone)]
pub struct Candidate {
    base_item : ItemDict,
    /// jaccard distance computed by hnsw
    distance : f32,
    ///
    likelyhood : f64,    
}


impl Candidate {
    pub fn new(item: ItemDict, distance:f32) -> Self {
        Candidate{ base_item : item, distance : distance, likelyhood : 0.}
    }

} // end of impl block Candidate

//==================================

/// A match between 2 sequences.
/// TODO devise an order
#[derive(Clone)]
pub struct SequenceMatch {
    request_item : ItemDict,
    ///
    candidates : Vec<Candidate>
}  // end of MatchSequence



impl SequenceMatch {
    pub fn new(request_item : ItemDict, candidates : Vec<Candidate>) -> Self {
        SequenceMatch{ request_item, candidates}
    }

}  // end of SequenceMatch


//====================================================================


pub struct GenomeMatch {
    /// filename containing request genome we want to identify
    request_file : String,
    /// target file in which some sequences were found 
    target_file : String,
    /// a sequence of matches between the request and a genome in the database.
    seq_matches : Vec<SequenceMatch>
}



/// A type providing keyed access to matches between requests and database genomes.
type Matches<'a> = HashMap<(String, String), GenomeMatch >;


/// This structure quantifies a match
#[derive(Clone)]
pub struct Matcher<'a> {
    kmer_size : usize,
    /// Dictionary of Database
    seqdict : &'a SeqDict,
    /// total number of bases of database.
    database_size : usize,
    ///
    seq_matches : Vec<SequenceMatch>,
}  // end of Matcher



impl <'a> Matcher<'a> {
    pub fn new(kmer_size : usize , seqdict : &'a SeqDict) -> Self {
        let database_size = seqdict.get_total_length();
        let seq_matches = Vec::with_capacity(1000);
        Matcher{kmer_size, seqdict, database_size, seq_matches}
    }

    pub fn insert_sequence_match(&mut self, request_item: ItemDict, candidates : Vec<Candidate>) {
        let new_match = SequenceMatch::new(request_item, candidates);
        self.seq_matches.push(new_match);
    }  // end of insert_sequence_match

    pub fn get_nb_sequence_match(&self) -> usize {
        self.seq_matches.len()
    }
    
} // end of impl block for Matcher