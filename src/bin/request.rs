//! Module request
//! try to match fasta sequence with repsect to database
//! 
//! request -- database [-d] dirname --query [-q] fastafileList --nbsearch [-n] nbanswers
//! 
//! - database is the name of directory ontaining hnsw dump files and seqdict dump
//! - fastaFileList is a file containing list of fasta file containing sequence to search for
//! - 
//



use clap::{App, Arg};

// for logging (debug mostly, switched at compile time in cargo.toml)
use env_logger::{Builder};

//
use std::time::{SystemTime};
use cpu_time::ProcessTime;

use archaea::idsketch::{SeqDict,Id};

fn main() {

}