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


use std::collections::HashMap;

use std::path::{PathBuf};
use std::fs::{OpenOptions};
use std::io::{Write,BufWriter};

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
    _likelyhood : f64,    
}


impl  SequenceMatch {
    pub fn new(base_item : ItemDict, distance:f32) -> Self {
        SequenceMatch{base_item : base_item, distance : distance, _likelyhood : 0.}
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
    // as far as we use split sequence coming from our split we do not need length as it is the same for all sequences
    fn compute_merit_wl(&self, threshold : f32) -> f64 {
        let mut merit = 1f64;
        for candidate in &self.candidates {
            if candidate.distance < threshold {
                merit *= candidate.distance as f64;
            }
        }
        return merit;
    }
}  // end of MatchList

//====================================================================

type GenomeMatch = HashMap<TargetGenome, MatchList>;

struct GenomeMatchAnalyzer {
} // end of GenomeMatchAnalyzer


impl <'a> GenomeMatchAnalyzer{

    pub fn new() -> Self {
        GenomeMatchAnalyzer{}
    }

    /// 
    /// returns vecors of (Id,similrity measure. The lesser the better).
    fn analyze(&mut self, targets : &GenomeMatch, threshold: f32) -> Vec<(String, f32)> {
        //
        let mut sorted = Vec::<(String,f32)>::new();
        // iterate through MatchList and call compute_merit
        for (g, m) in targets.iter() {
            let merit = m.compute_merit_wl(threshold); // the lower the better
            sorted.push((g.to_string(), merit as f32));
        }
        // sort in increasing order so first is better for us!
        sorted.sort_unstable_by(|a,b| a.1.partial_cmp(&b.1).unwrap());
        //
        sorted
    }

} // end of impl GenomeMatchAnalyzer

//===============================================================

// compute each genome length. Useful only in sequence mode not in block genome 
fn compute_genome_length(sequences : &SeqDict) -> HashMap<String, usize> {
    let mut hashed_length = HashMap::<String, usize>::new();
    //
    for seq in &sequences.0 {
        let entry = hashed_length.entry(seq.get_id().get_path().to_string()).or_insert(0);
        *entry += seq.get_len();
    }
    //
    return hashed_length;
}


/// This structure gather all sequence matches collected in function sketch_and_request_dir_compressedkmer.
/// It must dispatch matches to GenomeMatch quantifies a match
#[derive(Clone)]
pub struct Matcher {
    /// kmer size in ketching
    _kmer_size : usize,
    /// size of sketch
    _sketch_size : usize,
    /// genomes length. Maps a file path to a length. Useful only in split mode.
#[allow(unused)]
    genome_length : HashMap<String, usize>,
    /// total number of bases of database.
    _database_size : usize,
    /// for each request genome, we maintain a 
    seq_matches : HashMap<RequestGenome,  HashMap<TargetGenome, MatchList> >,
    ///
    nb_sequence_match : usize,
}  // end of Matcher



impl Matcher{
    pub fn new(_kmer_size : usize , _sketch_size: usize, seqdict : &SeqDict) -> Self {
        let _database_size = seqdict.get_total_length();
        let genome_length = compute_genome_length(seqdict);
        let seq_matches = HashMap::<RequestGenome, HashMap<TargetGenome, MatchList> >::new();
        Matcher{_kmer_size, _sketch_size, genome_length, _database_size, seq_matches, nb_sequence_match : 0}
    }



    pub fn insert_sequence_match(&mut self, req_item : ItemDict, new_matches : Vec<SequenceMatch>) {
        //
        // find request genome, if new we have to add an entry into seq_matches
        //                      if not get the hashmap of its TargetGenome and add entry.
        //
        let req_genome_path = req_item.get_id().get_path();
        let genome_for_insertion;
        if self.seq_matches.get_mut(req_genome_path).is_none() {
            // if request is new, insert it in request hashmap
            log::trace!("insert_sequence_match , creating entry for request genome {}", req_genome_path);
            self.seq_matches.insert(req_genome_path.clone(), HashMap::<TargetGenome, MatchList>::new());
        }
        genome_for_insertion = self.seq_matches.get_mut(req_genome_path).unwrap();
        // now we can add new_matches in genome_for_insertion
        for new_match in &new_matches {
            let candidate_item = new_match.get_item();
            let candidate_path = new_match.get_path();
            if  genome_for_insertion.get(candidate_path).is_none() {
             // do we have already candidate_genome in candidate list ?, if not insert it in target list
             log::trace!("insert_sequence_match , creating entry for target genome {}", candidate_path);
             genome_for_insertion.insert(candidate_path.to_string(), MatchList::new(candidate_item.clone()));
            }
            let match_list = genome_for_insertion.get_mut(candidate_path).unwrap();
            // insert candidate
            match_list.insert(new_match);
        } // end of for on new matches
        //
        self.nb_sequence_match += new_matches.len();
    }  // end of insert_sequence_match



    /// return mut reference to Target genome hashmap  matched for a given request genome, None if request is seen for the first time. 
    pub fn get_genome_match_mut(&mut self, request_path : &RequestGenome) -> Option<&mut HashMap<TargetGenome, MatchList> > {
        self.seq_matches.get_mut(request_path)
    }

    /// return total number of sequence matched. (useful only in sequence matching mode, not in one block genome match)
    pub fn get_nb_sequence_match(&self) -> usize {
        self.nb_sequence_match
    }
    
    /// This function must order all target genome match according a likelyhood/merit function
    /// The merit function is the sum on matched sequences for each target genome of (1-distance) * sequence length / total genome length.
    /// So it the fraction of length matched. The larger the better.
    pub fn analyze(&mut self) -> Result<(), String> {
        //
        let threshold = 0.99; // TODO ...
        let outname = "archea.matches";
        let outpath = PathBuf::from(outname.clone());
        let outfile = OpenOptions::new().write(true).create(true).truncate(true).open(&outpath);
        if outfile.is_err() {
            log::error!("Matcher analyze : dump could not open file {:?}", outpath.as_os_str());
            println!("Matcher analyze: could not open file {:?}", outpath.as_os_str());
            return Err("Opening of Matcher output failed").unwrap();
        }
        let mut outfile = BufWriter::new(outfile.unwrap());
        log::info!("dumping sorted match in : {}, threshold dist : {} ", outname, threshold);
        let mut match_analyzer = GenomeMatchAnalyzer::new();
        // iterate on request genome request
        for (genome, candidates) in self.seq_matches.iter_mut() {
            let sorted_match = match_analyzer.analyze(candidates, threshold);
            // print
            write!(outfile, "\n\n request genome : {}", genome).unwrap();
            let max_out = 5.min(sorted_match.len());
            for i in 0..max_out  {
                write!(outfile, "\n\t matched genome {}  merit : {:.3E}", sorted_match[i].0, sorted_match[i].1).unwrap();
            }
        }
        Ok(())
    } // end of analyze



} // end of impl block for Matcher


//===========================================================================

