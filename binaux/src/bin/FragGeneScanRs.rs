//! FragGeneScanRs executable
#![allow(non_snake_case)]

use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::Mutex;
use clap::{Arg, Command};
use seq_io::fasta;
use rayon::iter::{ParallelBridge, ParallelIterator};
use rayon::ThreadPoolBuilder;
use frag_gene_scan_rs::dna::{count_cg_content, Nuc};
use frag_gene_scan_rs::hmm;
use frag_gene_scan_rs::viterbi::viterbi;
use anyhow::Result;
use env_logger::Builder;

pub fn init_log() -> u64 {
    Builder::from_default_env().init();
    println!("\n ************** initializing logger *****************\n");
    return 1;
}

fn main() -> Result<()> {
    let _ = init_log();
    let matches = Command::new("FragGeneScanRs")
        .version("0.3.2")
        .about("Scalable high-throughput short-read open reading frame prediction.")
        .arg(Arg::new("seq-file")
             .short('s')
             .long("seq-file-name")
             .value_name("seq_file_name")
             .default_value("stdin")
             .help("Sequence file name including the full path. Using 'stdin' (or not supplying this argument) reads from standard input.")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("output-prefix")
             .short('o')
             .long("output-prefix")
             .value_name("output_prefix")
             .help("Output metadata (.out and .gff), proteins (.faa) and genes (.ffn) to files with this prefix. Use 'stdout' to write the predicted proteins to standard output.")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("complete")
             .short('w')
             .long("complete")
             .value_name("complete")
             .default_value("0")
             .help("The input sequence has complete genomic sequences; not short sequence reads.")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("formatted")
             .short('f')
             .long("formatted")
             .help("Format the DNA output.")
             .action(clap::ArgAction::SetTrue))
        .arg(Arg::new("train-file")
             .short('t')
             .long("training-file")
             .value_name("train_file_name")
             .required(true)
             .help("File name that contains model parameters; this file should be in the -r directory or one of the following:
             [complete] for complete genomic sequences or short sequence reads without sequencing error\n
             [sanger_5] for Sanger sequencing reads with about 0.5% error rate
             [sanger_10] for Sanger sequencing reads with about 1% error rate
             [454_5] for 454 pyrosequencing reads with about 0.5% error rate
             [454_10] for 454 pyrosequencing reads with about 1% error rate
             [454_30] for 454 pyrosequencing reads with about 3% error rate
             [illumina_1] for Illumina sequencing reads with about 0.1% error rate
             [illumina_5] for Illumina sequencing reads with about 0.5% error rate
             [illumina_10] for Illumina sequencing reads with about 1% error rate")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("train-file-dir")
             .short('r')
             .long("train-file-dir")
             .value_name("train_file_dir")
             .help("Full path of the directory containing the training model files.")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("thread-num")
             .short('p')
             .long("thread-num")
             .value_name("thread_num")
             .default_value("1")
             .help("The number of threads used by FragGeneScan++.")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("meta-file")
             .short('m')
             .long("meta-file")
             .value_name("meta_file")
             .help("Output metadata to this file (supersedes -o).")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("gff-file")
             .short('g')
             .long("gff-file")
             .value_name("gff_file")
             .help("Output metadata to this gff formatted file (supersedes -o).")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("aa-file")
             .short('a')
             .long("aa-file")
             .value_name("aa_file")
             .help("Output predicted proteins to this file (supersedes -o).")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("nucleotide-file")
             .short('n')
             .long("nucleotide-file")
             .value_name("nucleotide_file")
             .help("Output predicted genes to this file (supersedes -o).")
             .action(clap::ArgAction::Set))
        .arg(Arg::new("unordered")
             .short('u')
             .long("unordered")
             .help("Do not preserve record order in output (faster).")
             .action(clap::ArgAction::SetTrue))
        .get_matches();

    let (global, locals) = hmm::get_train_from_file(
        PathBuf::from(matches.get_one::<String>("train-file-dir").unwrap_or(&"train".to_string())),
        PathBuf::from(matches.get_one::<String>("train-file").unwrap()),
    )?;

    let input_file = matches.get_one::<String>("seq-file").map(String::as_str).unwrap_or("stdin");
    let input_seqs: Box<dyn Read + Send> = match input_file {
        "stdin" => Box::new(io::stdin()),
        filename => Box::new(File::open(filename)?),
    };

    let output_prefix = matches.get_one::<String>("output-prefix").map(String::as_str);

    let mut aastream: Option<Box<dyn Write + Send>> = match (
        matches.get_one::<String>("aa-file").map(String::as_str),
        output_prefix
    ) {
        (Some(filename), _) => Some(Box::new(File::create(filename)?)),
        (None, Some("stdout")) => Some(Box::new(io::stdout())),
        (None, Some(prefix)) => Some(Box::new(File::create(format!("{}.faa", prefix))?)),
        (None, None) => None,
    };

    let metastream: Option<Box<dyn Write + Send>> = match (
        matches.get_one::<String>("meta-file").map(String::as_str),
        output_prefix
    ) {
        (Some(filename), _) => Some(Box::new(File::create(filename)?)),
        (None, Some("stdout")) => None,
        (None, Some(prefix)) => Some(Box::new(File::create(format!("{}.out", prefix))?)),
        (None, None) => None,
    };

    let gffstream: Option<Box<dyn Write + Send>> = match (
        matches.get_one::<String>("gff-file").map(String::as_str),
        output_prefix
    ) {
        (Some(filename), _) => Some(Box::new(File::create(filename)?)),
        (None, Some("stdout")) => None,
        (None, Some(prefix)) => Some(Box::new(File::create(format!("{}.gff", prefix))?)),
        (None, None) => None,
    };

    let dnastream: Option<Box<dyn Write + Send>> = match (
        matches.get_one::<String>("nucleotide-file").map(String::as_str),
        output_prefix
    ) {
        (Some(filename), _) => Some(Box::new(File::create(filename)?)),
        (None, Some("stdout")) => None,
        (None, Some(prefix)) => Some(Box::new(File::create(format!("{}.ffn", prefix))?)),
        (None, None) => None,
    };

    if aastream.is_none() && metastream.is_none() && gffstream.is_none() && dnastream.is_none() {
        aastream = Some(Box::new(io::stdout()));
    }

    if matches.contains_id("unordered") {
        run(
            global,
            locals,
            input_seqs,
            aastream.map(UnbufferingBuffer::new),
            metastream.map(UnbufferingBuffer::new),
            gffstream.map(UnbufferingBuffer::new),
            dnastream.map(UnbufferingBuffer::new),
            matches.get_one::<String>("complete").map(|s| s == "1").unwrap_or(false),
            matches.contains_id("formatted"),
            matches.get_one::<String>("thread-num").and_then(|s| s.parse::<usize>().ok()).unwrap_or(1),
        )?;
    } else {
        run(
            global,
            locals,
            input_seqs,
            aastream.map(SortingBuffer::new),
            metastream.map(SortingBuffer::new),
            gffstream.map(SortingBuffer::new),
            dnastream.map(SortingBuffer::new),
            matches.get_one::<String>("complete").map(|s| s == "1").unwrap_or(false),
            matches.contains_id("formatted"),
            matches.get_one::<String>("thread-num").and_then(|s| s.parse::<usize>().ok()).unwrap_or(1),
        )?;
    }

    Ok(())
}

fn run<R: Read + Send, W: WritingBuffer + Send>(
    global: Box<hmm::Global>,
    locals: Vec<hmm::Local>,
    inputseqs: R,
    aa_buffer: Option<W>,
    meta_buffer: Option<W>,
    gff_buffer: Option<W>,
    dna_buffer: Option<W>,
    whole_genome: bool,
    formatted: bool,
    thread_num: usize,
) -> Result<()> {
    ThreadPoolBuilder::new()
        .num_threads(thread_num)
        .build_global()?;

    let meta_buffer = meta_buffer.map(Mutex::new);
    let gff_buffer = gff_buffer.map(Mutex::new);
    let dna_buffer = dna_buffer.map(Mutex::new);
    let aa_buffer = aa_buffer.map(Mutex::new);

    Chunked::new(100, fasta::Reader::new(inputseqs).into_records())
        .enumerate()
        .par_bridge()
        .map(|(index, recordvec)| {
            let mut metabuf = Vec::new();
            let mut gffbuf = Vec::new();
            let mut dnabuf = Vec::new();
            let mut aabuf = Vec::new();
            for record in recordvec {
                let fasta::OwnedRecord { mut head, seq } = record?;
                head = head.into_iter().take_while(u8::is_ascii_graphic).collect();
                let nseq: Vec<Nuc> = seq.into_iter().map(Nuc::from).collect();
                let read_prediction = viterbi(
                    &global,
                    &locals[count_cg_content(&nseq)],
                    head,
                    nseq,
                    whole_genome,
                );
                if meta_buffer.is_some() {
                    read_prediction.meta(&mut metabuf)?;
                }
                if gff_buffer.is_some() {
                    read_prediction.gff(&mut gffbuf)?;
                }
                if dna_buffer.is_some() {
                    read_prediction.dna(&mut dnabuf, formatted)?;
                }
                if aa_buffer.is_some() {
                    read_prediction.protein(&mut aabuf, whole_genome)?;
                }
            }
            if let Some(buffer) = &meta_buffer {
                buffer.lock().unwrap().add(index, metabuf)?;
            }
            if let Some(buffer) = &gff_buffer {
                buffer.lock().unwrap().add(index, gffbuf)?;
            }
            if let Some(buffer) = &dna_buffer {
                buffer.lock().unwrap().add(index, dnabuf)?;
            }
            if let Some(buffer) = &aa_buffer {
                buffer.lock().unwrap().add(index, aabuf)?;
            }
            Ok(())
        })
        .collect()
}

struct Chunked<I: Iterator> {
    size: usize,
    iterator: I,
}

impl<I: Iterator> Chunked<I> {
    fn new(size: usize, iterator: I) -> Self {
        Chunked { size, iterator }
    }
}

impl<I: Iterator> Iterator for Chunked<I> {
    type Item = Vec<I::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut items = Vec::with_capacity(self.size);
        for _ in 0..self.size {
            if let Some(item) = self.iterator.next() {
                items.push(item);
            } else {
                break;
            }
        }
        if items.is_empty() {
            None
        } else {
            Some(items)
        }
    }
}

trait WritingBuffer {
    fn add(&mut self, index: usize, item: Vec<u8>) -> Result<()>;
}

struct SortingBuffer<W: Write + Send> {
    next: usize,
    queue: VecDeque<Option<Vec<u8>>>,
    stream: W,
}

impl<W: Write + Send> SortingBuffer<W> {
    fn new(stream: W) -> Self {
        SortingBuffer {
            next: 0,
            queue: VecDeque::new(),
            stream: stream,
        }
    }
}

impl<W: Write + Send> WritingBuffer for SortingBuffer<W> {
    fn add(&mut self, index: usize, item: Vec<u8>) -> Result<()> {
        while self.next + self.queue.len() <= index {
            self.queue.push_back(None);
        }
        self.queue[index - self.next] = Some(item);

        while self.queue.front().map(Option::is_some).unwrap_or(false) {
            let item = self.queue.pop_front().unwrap().unwrap();
            self.next += 1;
            self.stream.write_all(&item)?;
        }
        Ok(())
    }
}

struct UnbufferingBuffer<W: Write + Send> {
    stream: W,
}

impl<W: Write + Send> UnbufferingBuffer<W> {
    fn new(stream: W) -> Self {
        UnbufferingBuffer { stream }
    }
}

impl<W: Write + Send> WritingBuffer for UnbufferingBuffer<W> {
    fn add(&mut self, _: usize, item: Vec<u8>) -> Result<()> {
        self.stream.write_all(&item)?;
        Ok(())
    }
}
