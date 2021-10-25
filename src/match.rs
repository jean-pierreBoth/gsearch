//!  A module to quantify match between sequence
//! For each file/genome we maintain a list of matches with other files/genomes
//! A match between sequences is a couple of sequence one from request genome, one from the database
//! for which we found a distance below some significant threshold.
//! 
//! A match between genomes is a list of matches between couple of sequences, one of the request genome, 
//! the other from the database genomes.


/// A match between 2 sequences.
/// TODO devise an order
pub SequenceMatch {
    /// fasta_id of request
    req_fasta_id : &String,
    /// length of request sequence
    req_len : usize,
    /// fasta_id of potential target
    target_fasta_id : &String,
    /// length of potential taret
    target_len : usize,
    /// jaccard distance computed by hnsw
    distance : f32,
    //
}  // end of MatchSequence



pub struct GenomeMatch {
    /// filename containing request genome we want to identify
    request_file : String,
    /// target file in which some sequences were found 
    target_file : String,
    /// a sequence of mathes between the request and a genome in the database.
    seq_matches : Vec<SequenceMatch>
}


/// A type providing keyed access to matches between requests and database genomes.
type Matches = HashMap<(String, String), GenomeMatch>