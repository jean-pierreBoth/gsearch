![alt text](https://github.com/jean-pierreBoth/archaea/blob/master/GSearch-logo.jpg?raw=true)

# A rust classifier based on probminhash and HNSW for microbial genomes

ARCHAEA stands for: <u>A</u> <u>R</u>ust <u>C</u>lassifier base on <u>H</u>ierarchical N<u>a</u>vigable SW graphs, <u>e</u>t.<u>a</u>l.** Later on, we rename it to GSearch, stands of Genomic Search.

This package (**currently in development**) compute probminhash signature of  bacteria and archaea (or virus and fungi) genomes and stores the id of bacteria and probminhash signature in a Hnsw structure for searching of new request genomes.

This package is developped by Jean-Pierre Both (https://github.com/jean-pierreBoth) for the software part and Jianshu Zhao (https://github.com/jianshu93) for the genomics part.

## Sketching of genomes/tohnsw

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
tohnsw -d db_dir_nt -s 12000 -k 16 --ef 1600 -n 128
tohnsw -d db_dir_aa -s 12000 -k 7 --ef 1600 -n 128 --aa
### request neighbours for each genomes (fna, fasta, faa et.al. are supported) in query_dir_nt or aa using pre-built database:
wget http://enve-omics.ce.gatech.edu/data/public_gsearch/GTDB_r207_hnsw_graph.tar.gz
tar xzvf ./GTDB_r207_hnsw_graph.tar.gz
cd ./GTDB_r207_hnsw_graph/nucl
### request neighbors for nt genomes
request -b ./ -d query_dir_nt -n 50
### request neighbors for aa genomes (predicated by Prodigal or FragGeneScanRs)
cd ./GTDB_r207_hnsw_graph/prot
request -b ./ -d query_dir_aa -n 50 --aa
### request neighbors for aa universal gene (extracted by hmmer according to hmm files provided)
cd ./GTDB_r207_hnsw_graph/universal
request -b ./ -d query_dir_universal_aa -n 50 --aa
```


## Dependencies, features and Installation

### features


* hnsw_rs relies on the crate simdeez to accelerate distance computation. On intel you can build hnsw_rs with the feature simdeez_f

* annembed relies on openblas so you must choose between  the features "annembed_openblas-static" , "annembed_openblas-system" or "annembed_intel-mkl". You may need to install gcc, gfortran and make.

* kmerutils provides a feature "withzmq". This feature can be used to store compressed qualities on a server and run requests. It is not necessary in this crate.

### Simple case for install:

**Pre-built binaries** are available on release page (https://github.com/jean-pierreBoth/archaea/releases/tag/v1.0) for major platforms. If you wan to install/compiling by yourself:

```bash
###A simple installation, with annembed enabled would be:
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install archaea --features="annembed_intel-mkl"

###on MacOS, which requires dynamic library link:
cargo build --release --features="annembed_openblas-system" 

##or on intel using openblas instead of intel-mkl:  
cargo build --release --features="annembed_openblas-system" --features="hnsw_rs/simdeez_f"

###Then install FragGeneScanRs:
cargo install --git https://gitlab.com/Jianshu_Zhao/fraggenescanrs
```

Alternatively it is possible to modify the features section in  Cargo.toml. Just fill in the default you want.

### Some hints in case of problem (including installing/compiling on ARM CPUs) are given [here](./installpb.md)

### Pre-built databases

We provide pre-built genome/proteome database graph file for bacteria/archaea, virus and fungi. Proteome database are based on genes for each genome, predicted by FragGeneScanRs (https://gitlab.com/Jianshu_Zhao/fraggenescanrs) for bacteria/archaea/virus and GeneMark-ES version 2 (http://exon.gatech.edu/GeneMark/license_download.cgi) for fungi.  

- Bacteria/archaea genomes are the newest version of GTDB database (https://gtdb.ecogenomic.org), which defines a bacterial speces at 95% ANI. Note that GSearch can also run for even higher resolution species database such as 99% ANI.
- Virus data base are based on the JGI IMG/VR database newest version (https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html), which also define a virus OTU (vOTU) at 95% ANI.  
- Fungi database are based on the entire RefSeq fungal genomes (retrived via the MycoCosm website), we dereplicated and define a fungal speices at 99.5% ANI. 
- All three pre-built database can be available here:http://enve-omics.ce.gatech.edu/data/gsearch 