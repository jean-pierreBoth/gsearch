//! A module to quantify match between sequence.  
//! 
//! For each file/genome we maintain a list of matches with other files/genomes
//! A match between sequences is a couple of sequence one from request genome, one from the database
//! for which we found a distance below some significant threshold.
//! 
//! A match between genomes is a list of matches between couple of sequences, one of the request genome, 
//! the other from the database genomes.  
//! 
//! Knowing sequence length of request and target in database, kmer size and sketch size we must get
//! a probability of match getting ideas from :
//! *How independent are the appearances of n-mers in different genomes Fufanov and al. BioInformatics 2004*

#![allow(unused)]

use std::collections::HashMap;

use crate::utils::*;


type RequestGenome = String;
type TargetGenome = String;


#[derive(Clone)]
pub struct SequenceMatch {
    /// sequence id in database.
    base_item : ItemDict,
    /// jaccard distance computed by hnsw
    distance : f32,
    ///
    likelyhood : f64,    
}


impl  SequenceMatch {
    pub fn new(base_item : ItemDict, distance:f32) -> Self {
        SequenceMatch{base_item : base_item, distance : distance, likelyhood : 0.}
    }

    /// return genome path of candidate
    pub fn get_path(&self) -> &String {
        self.base_item.get_id().get_path()
    }

    pub fn get_item(&self) -> &ItemDict {
        &self.base_item
    }

} // end of impl block Candidate

//==================================

/// A candidate genome stores sequence matches involving it 
/// So all candidates in this structure belong to th same candidate genome.
/// TODO devise an order
#[derive(Clone)]
pub struct MatchList {
    /// 
    base_item : ItemDict,
    ///
    candidates : Vec<SequenceMatch>
}  // end of MatchSequence



impl MatchList  {
    pub fn new(base_item : ItemDict) -> Self {

        MatchList { base_item, candidates : Vec::<SequenceMatch>::new()}
    }

    pub fn get_item(&self) -> &ItemDict {
        &self.base_item
    }

    pub fn get_path(&self) -> &String {
        &self.base_item.get_id().get_path()
    }

    fn insert(&mut self, seq_match : &SequenceMatch) {
        self.candidates.push(seq_match.clone());
    }
}  // end of MatchList


//====================================================================



type GenomeMatch =  HashMap<TargetGenome, MatchList>;

/// This structure gther all sequence matches collected in function sketch_and_request_dir_compressedkmer.
/// It must dispatch matches to GenomeMatch quantifies a match
#[derive(Clone)]
pub struct Matcher<'a> {
    /// kmer size in ketching
    kmer_size : usize,
    /// size of sketch
    sketch_size : usize,
    /// Dictionary of Database
    seqdict : &'a SeqDict,
    /// total number of bases of database.
    database_size : usize,
    /// for each request genome, we maintain a 
    seq_matches : HashMap<RequestGenome,  HashMap<TargetGenome, MatchList> >,
    ///
    nb_sequence_match : usize,
}  // end of Matcher



impl <'a> Matcher<'a> {
    pub fn new(kmer_size : usize , sketch_size: usize, seqdict : &'a SeqDict) -> Self {
        let database_size = seqdict.get_total_length();
        let seq_matches = HashMap::<RequestGenome, HashMap<TargetGenome, MatchList> >::new();
        Matcher{kmer_size, sketch_size, seqdict, database_size, seq_matches, nb_sequence_match : 0}
    }

    pub fn insert_sequence_match(&mut self, req_item : ItemDict, new_matches : Vec<SequenceMatch>) {
        //
        // find request genome, if new we have to add an entry into seq_matches
        //                      if not get the hashmap of its TargetGenome and add entry.
        //
        let req_genome_path = req_item.get_id().get_path();
        let genome_match = self.seq_matches.get_mut(req_genome_path);
        let mut genome_for_insertion;
        if self.seq_matches.get_mut(req_genome_path).is_none() {
            // if request is new, insert it in request hashmap
            log::info!("insert_sequence_match , creating entry for request genome {}", req_genome_path);
            self.seq_matches.insert(req_genome_path.clone(), HashMap::<TargetGenome, MatchList>::new());
        }
        genome_for_insertion = self.seq_matches.get_mut(req_genome_path).unwrap();
        // now we can add new_matches in genome_for_insertion
        for new_match in &new_matches {
            let candidate_item = new_match.get_item();
            let candidate_path = new_match.get_path();
            if  genome_for_insertion.get(candidate_path).is_none() {
             // do we have already candidate_genome in candidate list ?, if not insert it in target list
             log::info!("insert_sequence_match , creating entry for target genome {}", req_genome_path);
             genome_for_insertion.insert(candidate_path.to_string(), MatchList::new(candidate_item.clone()));
            }
            let match_list = genome_for_insertion.get_mut(candidate_path).unwrap();
            // insert candidate
            match_list.insert(new_match);
        } // end of for on new matches
        //
        self.nb_sequence_match += new_matches.len();
    }  // end of insert_sequence_match

    pub fn get_genome_match(&self, request : &RequestGenome) -> Option<&HashMap<TargetGenome, MatchList> > {
        None
    }

    pub fn get_nb_sequence_match(&self) -> usize {
        self.nb_sequence_match
    }
    
    /// This function must order all target genome match according a likelyhood/merit function
    /// The merit function is the sum on matched sequences for each target genome of (1-distance) * sequence length / total genome length. 
    pub fn analyze(&self) {

    }
} // end of impl block for Matcher


//===========================================================================

