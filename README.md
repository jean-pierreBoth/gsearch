![alt text](https://github.com/jean-pierreBoth/archaea/blob/master/GSearch-logo.jpg?raw=true)

# A rust classifier based on probminhash and HNSW for prokaryotic genomes

ARCHAEA stands for: <u>A</u> <u>R</u>ust <u>C</u>lassifier base on <u>H</u>ierarchical N<u>a</u>vigable SW graphs, <u>e</u>t.<u>a</u>l.**

This package (**currently in development**) compute probminhash signature of  bacteria and archaea genomes and stores the id of bacteria and probminhash signature in a Hnsw structure for searching of new request genomes.

This package is developped by Jean-Pierre Both (https://github.com/jean-pierreBoth) for the software part and Jianshu Zhao (https://github.com/jianshu93) for the genomics part.

## Sketching of genomes

The sketching and database is done by the module ***tohnsw***.

The sketching of reference genomes can take some time (one or 2 hours for 50000 bacterial genomes of NCBI for parameters giving a correct quality of sketching). Result is stored in 2 structures:
- A Hnsw structure storing rank of data processed and corresponding sketches.
- A Dictionary associating each rank to a fasta id and fasta filename.

The Hnsw structure is dumped *in hnswdump.hnsw.graph* and  *hnswdump.hnsw.data*
The Dictionary is dumped in a json file *seqdict.json*
## Requests

For requests  the module ***request*** is being used. It reloads the dumped files, hnsw and seqdict related
takes a list of fasta files containing requests and for each fasta file dumps the asked number of nearest neighbours.

## Usage

```bash
### build database given genome file directory, fna.gz was expected. L for nt and .faa or .faa.gz for --aa. Limit for k is 32 (15 not work due to compression), for s is 65535 (u16) and for n is 255 (u8)
tohnsw -d db_dir_nt -s 12000 -k 21 --ef 1600 -n 128
tohnsw -d db_dir_aa -s 24000 -k 7 --ef 1600 -n 128 --aa
### request neighbours for each genomes in query_dir given pre-built database path
request -b ./ -d query_dir_nt -n 50
request -b ./ -d query_dir_aa -n 50 --aa
```


## Dependencies, features and Installation

### features


* hnsw_rs relies on the crate simdeez to accelerate distance computation. On intel you can build hnsw_rs with the feature simdeez_f

* annembed relies on openblas so you must choose between  the features "annembed_openblas-static" , "annembed_openblas-system" or "annembed_intel-mkl". You may need to install gcc, gfortran and make.

* kmerutils provides a feature "withzmq". This feature can be used to store compressed qualities on a server and run requests. It is not necessary in this crate.

### simple case

    A simple installation, with annembed enabled would be:
    cargo install 
    cargo build --release --features="annembed_openblas-system" 
    or on intel :  
    cargo build --release --features="annembed_openblas-system" --features="hnsw_rs/simdeez_f"

Alternatively it is possible to modify the features section in  Cargo.toml. Just fill in the default you want.

### Some hints in case of problem are given [here](./installpb.md)
### Homology search




The last step involves a homology search using hmmer, which can be directly installed using conda or brew. If you are using apple M1 ARM/aarch64 structure, you can have a native support of hmmer folloing the steps:

```bash
### download h3-heno branch of hmmer here (do not git clone but download zip):

https://github.com/EddyRivasLab/hmmer/tree/h3-arm

## go into the donwloaded directory and download Easel develop branch here (do not git clone but download zip) :
cd h3-arm
https://github.com/EddyRivasLab/easel/tree/develop

## compile, or you can download binaries from here: https://github.com/jianshu93/hmmer-h3-arm
autoconf
./configure
make -j 8
sudo make install
hmmsearch -h
```


### Pre-built databases

We provide pre-built genome/proteome database graph file for bacteria/archaea, virus and fungi. Proteome database are based on genes for each genome, predicted by FragGeneScanRs (https://gitlab.com/Jianshu_Zhao/fraggenescanrs) for bacteria/archaea/virus and GeneMark-ES version 2 (http://exon.gatech.edu/GeneMark/license_download.cgi) for fungi.  

- Bacteria/archaea genomes are the newest version of GTDB database (https://gtdb.ecogenomic.org), which defines a bacterial speces at 95% ANI. Note that ARCHAEA can also run for even higher resolution species database such as 99% ANI.

- Virus data base are based on the JGI IMG/VR database newest version (https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html), which also define a virus OTU (vOTU) at 95% ANI.  

- Fungi database are based on the entire RefSeq fungal genomes, we dereplicated and define a fungal speices at 99% ANI. All three pre-built database can be available here: