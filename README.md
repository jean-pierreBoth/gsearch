# Gsearch: A rust classifier based on MinHash metric and HNSW for microbial genomes

![Alt!](https://github.com/jean-pierreBoth/gsearch/blob/master/GSearch-logo.jpg?raw=true)

gsearch is the new name of the crate archaea.  It stands for **genomic search**.  

This package (**currently in development**) compute probminhash signature of  bacteria and archaea (or virus and fungi) genomes and stores the id of bacteria and probminhash signature in a Hnsw structure for searching of new request genomes.

This package is developped by Jean-Pierre Both [jpboth](https://github.com/jean-pierreBoth) for the software part and [Jianshu Zhao](https://github.com/jianshu93) for the genomics part. We also created a mirror of this repo on [GitLab](https://gitlab.com/Jianshu_Zhao/gsearch) and [Gitee](https://gitee.com/jianshuzhao/gsearch), just in case Github service is not available in some region, e.g, China.

## Sketching of genomes/tohnsw

The objective is to use the Jaccard index as an accurate proxy of mutation rate or Average Nucleitide Identity(ANI). To achieve this we use sketching.  
We generate kmers along sequences and sketch the kmer distribution encountered in a file. Then final sketch is stored in a Hnsw database See [hnsw](https://arxiv.org/abs/1603.09320).

The sketching and database is done by the subcommand ***tohnsw***.

The Jaccard index come in 2 flavours:  
    1. The probability Jaccard index that takes into account the Kmer multiplicity. It is defined by :
    $$J_{P(A,B)}=\sum_{d\in D} \frac{1}{\sum_{d'\in D} \max (\frac{\omega_{A}(d')}{\omega_{A}(d)},\frac{\omega_{B}(d')}{\omega_{B}(d)})}$$
    where $\omega_{A}(d)$ is the multiplicity of $d$ in A
    (see [Moulton-Jiang-arxiv](https://arxiv.org/abs/1809.04052)).
    In this case we use the probminhash algorithm as implemented in [probminhash](https://github.com/jean-pierreBoth/probminhash)  
    2. The unweighted (simple) Jaccard index defined by :
        $$Jaccard(A,B)=\frac{A \cap B}{A \cup B}$$
        In this case we use the SuperMinhash or the SetSketch (based on hyperloglog) method.  
    The estimated Jaccard index is used to build HNSW graph database, which is implemented in crate [hnswlib-rs](https://crates.io/crates/hnsw_rs).

The sketching of reference genomes can take some time (one or 2 hours for ~65,000 bacterial genomes of NCBI for parameters giving a correct quality of sketching). Result is stored in 2 structures:

- A Hnsw structure storing rank of data processed and corresponding sketches.
- A Dictionary associating each rank to a fasta id and fasta filename.

The Hnsw structure is dumped *in hnswdump.hnsw.graph* and  *hnswdump.hnsw.data*
The Dictionary is dumped in a json file *seqdict.json*

## Requests

For requests  the subcommand ***request*** is being used. It reloads the dumped files, hnsw and seqdict related
takes a list of fasta files containing requests and for each fasta file dumps the asked number of nearest neighbours.

## Simple case for install

**Pre-built binaries** will be available on release page [binaries](https://github.com/jean-pierreBoth/gsearch/releases/tag/v0.0.12) for major platforms (no need to install but just download and make it executable). We recommend you use the linux one (GSearch_Linux_x86-64_intel-mkl-static.zip
) for linux system in this release page for convenience because the only dependency is GCC (Recent Linux version does not allow static compiling of GCC libraries like libc.so.6). For macOS, we recommend the universal binary [mac-binaries](https://github.com/jean-pierreBoth/gsearch/releases/download/v0.0.12/GSearch_darwin_universal.zip) for any macOS platform (x86-64 or arm64).

## Or if you have conda installed

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gsearch/README.html)

```bash
    conda config --add channels bioconda
    conda install gsearch -c bioconda
```

Otherwise it is possible to install/compile by yourself (see install section)

```bash

### get the binary for linux (make sure you have recent Linux installed with GCC, e.g., Ubuntu 18.0.4 or above)

wget <https://github.com/jean-pierreBoth/gsearch/releases/download/v0.0.12/GSearch_Linux_x86-64_intel-mkl-static.zip>
unzip GSearch_Linux_x86-64_intel-mkl-static.zip

## get the binary for macOS

wget <https://github.com/jean-pierreBoth/gsearch/releases/download/v0.0.12/GSearch_darwin_universal.zip>
unzip GSearch_darwin_universal.zip

* **make it excutable (changed it accordingly on macOS)**

chmod a+x ./GSearch_Linux_x86-64_intel-mkl-static/*


### put it under your system/usr bin directory (/usr/local/bin/ as an example) where it can be called

mv ./GSearch_Linux_x86-64_intel-mkl-static/* /usr/local/bin/

### check install

tohnsw -h
request -h

### check install MacOS, you may need to change the system setup to allow external binary to run by type the following first and use your admin password

sudo spctl --master-disable

### and then

gsearchbin -h

```

## usage

We give here an example of utilization with prebuilt databases.

```bash
### download neighbours for each genomes (fna, fasta, faa et.al. are supported) in query_dir_nt or aa using pre-built database

wget <http://enve-omics.ce.gatech.edu/data/public_gsearch/GTDB_r207_hnsw_graph.tar.gz>
tar xzvf ./GTDB_r207_hnsw_graph.tar.gz

### get test data, we provide 2 genomes at nt, AA and universal gene level

wget <https://github.com/jean-pierreBoth/gsearch/releases/download/v0.0.12/test_data.tar.gz>
tar xzvf ./test_data.tar.gz

cd ./GTDB_r207_hnsw_graph/nucl

### request neighbors for nt genomes (here -n is how many neighbors you want to return for each of your query genome)

gsearchbin request -b ./ -r ../../test_data/query_dir_nt -n 50

### or request neighbors for aa genomes (predicted by Prodigal or FragGeneScanRs)

cd ./GTDB_r207_hnsw_graph/prot
gsearchbin request -b ./ -r ../../test_data/query_dir_aa -n 50 --aa

### or request neighbors for aa universal gene (extracted by hmmer according to hmm files from gtdb, we also provide one in release page)

cd ./GTDB_r207_hnsw_graph/universal
gsearchbin request -b ./ -r ../../test_data/query_dir_universal_aa -n 50 --aa

### Building database. database is huge in size, users are welcome to download gtdb database here: (<https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz>) and here (<https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/genomic_files_reps/gtdb_proteins_aa_reps_r207.tar.gz>)

### build database given genome file directory, fna.gz was expected. L for nt and .faa or .faa.gz for --aa. Limit for k is 32 (15 not work due to compression), for s is 65535 (u16) and for n is 255 (u8)

gsearchbin tohnsw -d db_dir_nt -s 12000 -k 16 --ef 1600 -n 128
gsearchbin tohnsw -d db_dir_aa -s 12000 -k 7 --ef 1600 -n 128 --aa

### When there are new genomes  after comparing with the current database (GTDB v207, e.g. ANI < 95% with any genome after searcing, corresponding to >0.875 ProbMinHash distance), those genomes can be added to the database

### must run in the existing database file folder

cd ./GTDB_r207_hnsw_graph/nucl

### old .graph,.data and all .json files will be updated to the new one. Then the new one can be used for requesting as an updated database

gsearchbin tohnsw -d db_dir_nt (new genomes directory) -s 12000 -k 16 --ef 1600 -n 128 --add

### or add at the amino acid level

cd ./GTDB_r207_hnsw_graph/prot
gsearchbin tohnsw -d db_dir_nt (new genomes directory in AA format predicted by prodigal/FragGeneScanRs) -s 12000 -k 16 --ef 1600 -n 128 --aa --add

```

## Output explanation

**gsearch.answers** is the default output file in your current directory.  
 For each of your genome in the query_dir, there will be requested N nearest genomes found and sorted by distance (smallest to largest).  
  If one genome in the query does not exist in the output file, meaning at this level (nt or aa), there is no such nearest genomes in the database (or distant away from the best hit in the database), you may then go to amino acid level or universal gene level.

## Dependencies, features and Installation

### Features

- hnsw_rs relies on the crate simdeez to accelerate distance computation. On intel you can build hnsw_rs with the feature simdeez_f

- annembed relies on openblas so you must choose between  the features "annembed_openblas-static" , "annembed_openblas-system" or "annembed_intel-mkl". You may need to install gcc, gfortran and make.
This can be done using the **--features** option as explained below, or by modifying the features section in  Cargo.toml. In that case just fill in the default you want.
- kmerutils provides a feature "withzmq". This feature can be used to store compressed qualities on a server and run requests. It is not necessary in this crate.

## Install

### First install Rust tools

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### gsearch installation and compilation from Crates.io (Recommended)

- simple installation, with annembed enabled would be with intel-mkl

```bash
    cargo install gsearch --features="annembed_intel-mkl"
```

or with a system installed openblas:

```bash
cargo install gsearch --features="annembed_openblas-system"
```

- On MacOS, which requires dynamic library link (you have to install openblas first and then xz, the MacOS/Darwin binary provided also requires this):
(note that openblas install lib path is different on M1 MACs).  
So you need to run:

```bash
    brew install openblas xz
    echo 'export LDFLAGS="-L/usr/local/opt/openblas/lib"' >> ~/.bash_profile
    echo 'export CPPFLAGS="-I/usr/local/opt/openblas/include"' >> ~/.bash_profile
    echo 'export PKG_CONFIG_PATH="/usr/local/opt/openblas/lib/pkgconfig"' >> ~/.bash_profile
    source ~/.bash_profile
    cargo install gsearch --features="annembed_openblas-system"
```

- Intel:  
  You can enable simd instruction with the feature hnsw_rs/simdeez_f.  
  Using openblas instead of intel-mkl you would run:

```bash
cargo install gsearch --features="annembed_openblas-system" --features="hnsw_rs/simdeez_f"
```

#### gsearch installation from the most recent version from github

- direct installation from github:
  
```bash
    cargo install gsearch --features="annembed_intel-mkl" --git https://github.com/jean-pierreBoth/gsearch
```

- download and compilation
  
```bash
git clone https://github.com/jean-pierreBoth/gsearch
cd gsearch
#### build
cargo build --release --features="annembed_openblas-static"
###on MacOS, which requires dynamic library link:
cargo build --release --features="annembed_openblas-system"
```

#### Documentation generation

Html documentation can be generated by running (example for someone using the "annembed_openblas-system" feature):

```bash
cargo doc --features="annembed_openblas-system" --no-deps --open
```

### Then install FragGeneScanRs

```bash
cargo install --git https://gitlab.com/Jianshu_Zhao/fraggenescanrs
```

### Some hints in case of problem (including installing/compiling on ARM64 CPUs) are given [here](./installpb.md)

## Pre-built databases

We provide pre-built genome/proteome database graph file for bacteria/archaea, virus and fungi. Proteome database are based on genes for each genome, predicted by FragGeneScanRs (<https://gitlab.com/Jianshu_Zhao/fraggenescanrs>) for bacteria/archaea/virus and GeneMark-ES version 2 (<http://exon.gatech.edu/GeneMark/license_download.cgi>) for fungi.  

- Bacteria/archaea genomes are the newest version of GTDB database (<https://gtdb.ecogenomic.org>), which defines a bacterial speces at 95% ANI. Note that GSearch can also run for even higher resolution species database such as 99% ANI.
- Virus data base are based on the JGI IMG/VR database newest version (<https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html>), which also define a virus OTU (vOTU) at 95% ANI.  
- Fungi database are based on the entire RefSeq fungal genomes (retrived via the MycoCosm website), we dereplicated and define a fungal speices at 99.5% ANI.
- All three pre-built databases are available here:<http://enve-omics.ce.gatech.edu/data/gsearch>

## References

1. Jianshu Zhao, Jean Pierre Both, Luis M. Rodriguez-R and Konstantinos T. Konstantinidis, 2022. GSearch: Ultra-Fast and Scalable Microbial Genome Search by combining Kmer Hashing with Hierarchical Navigable Small World Graphs. *bioRxiv* 2022:2022.2010.2021.513218. [biorxiv](https://www.biorxiv.org/content/10.1101/2022.10.21.513218v2).
