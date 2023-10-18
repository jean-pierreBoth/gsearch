# GSearch: A Rust Genomic Search Program based on Various MinHash-like Metrics and HNSW for Microbial Genomes

![Alt!](https://github.com/jean-pierreBoth/gsearch/blob/master/GSearch-logo.jpg?raw=true)

gsearch is the new name of the crate archaea.  It stands for **genomic search**.  

This package (**currently in development**) compute probminhash signature of  bacteria and archaea (or virus and fungi) genomes and stores the id of bacteria and probminhash signature in a Hnsw structure for searching of new request genomes.

This package is developped by Jean-Pierre Both [jpboth](https://github.com/jean-pierreBoth) for the software part and [Jianshu Zhao](https://github.com/jianshu93) for the genomics part. We also created a mirror of this repo on [GitLab](https://gitlab.com/Jianshu_Zhao/gsearch) and [Gitee](https://gitee.com/jianshuzhao/gsearch) (You need to log in first to see the content), just in case Github service is not available in some region, e.g, China.

## Sketching of genomes/tohnsw

The objective is to use the Jaccard index as an accurate proxy of mutation rate or Average Nucleitide Identity(ANI) or Average Amino Acide Identity (AAI) According to equation:
$$ANI=1+\frac{1}{k}log\frac{2*J}{1+J}$$

where J is Jaccard-like index (e.g. Jp from ProbMinHash or J from SuperMinHash, SetSketch or Densified MinHash, a class of locality sensitive hashing algorithms, suitable for nearest neighbor search,  see below) and k is k-mer size.
To achieve this we use sketching. We generate kmers along sequences and sketch the kmer distribution encountered in a file. Then final sketch is stored in a Hnsw database See [hnsw](https://arxiv.org/abs/1603.09320).

The sketching and database is done by the subcommand ***tohnsw***.

The Jaccard index come in 2 flavours:  
    1. The probability Jaccard index that takes into account the Kmer multiplicity. It is defined by :
    $$J_{P(A,B)}=\sum_{d\in D} \frac{1}{\sum_{d'\in D} \max (\frac{\omega_{A}(d')}{\omega_{A}(d)},\frac{\omega_{B}(d')}{\omega_{B}(d)})}$$
    where $\omega_{A}(d)$ is the multiplicity of $d$ in A
    (see [Moulton-Jiang-arxiv](https://arxiv.org/abs/1809.04052)).
    In this case for J_p we use the probminhash algorithm as implemented in [probminhash](https://github.com/jean-pierreBoth/probminhash)  
    2. The unweighted (simple) Jaccard index defined by :
        $$Jaccard(A,B)=\frac{A \cap B}{A \cup B}$$
        In this case for J we use the [SuperMinHash](https://arxiv.org/abs/1706.05698), [SetSketch](https://vldb.org/pvldb/vol14/p2244-ertl.pdf) (based on the locality sensitive hashing) or Densified [MinHash](http://proceedings.mlr.press/v70/shrivastava17a.html) method, also implemented in probminhash mentioned above.  
    The estimated Jaccard-like index is used to build HNSW graph database, which is implemented in crate [hnswlib-rs](https://crates.io/crates/hnsw_rs).

The sketching of reference genomes can take some time (less than one hours for ~65,000 bacterial genomes of GTDB for parameters giving a correct quality of sketching, or 3 to 4 hours for entire NCBI/RefSeq prokaryotic genomes. ~318K). Result is stored in 2 structures:

- A Hnsw structure storing rank of data processed and corresponding sketches.
- A Dictionary associating each rank to a fasta id and fasta filename.

The Hnsw structure is dumped *in hnswdump.hnsw.graph* and  *hnswdump.hnsw.data*
The Dictionary is dumped in a json file *seqdict.json*

## Adding genomes to existing database
For adding new genomes to existing database, the ***add*** subcommand is being used. It will automatically load sketching files, graph files and also paremeters used for building the graph and then use the same parameters to add new genomes to exisiting database genomes.

## Request

For requests  the subcommand ***request*** is being used. It reloads the dumped files, hnsw and seqdict related
takes a list of fasta files containing requests and for each fasta file dumps the asked number of nearest neighbours according to distance mentioned above. A tabular file will be saved to disk with 3 key columns: query genome path, database genome path (ranked by distance) and distance. The distance can be transformed into ANI or AAI according to the equation above. Check the scripts for analyzing output from request here: https://github.com/jianshu93/gsearch_analysis

## Ann
For UMAP-like algorithm to perform dimension reduction and then visuzlizing genome database, we run it after the tohnsw step (pre-built database) (see below useage ann section). See [annembed](https://github.com/jean-pierreBoth/annembed) crate for details. Then the output of this step can be visualized, for example for the GTDB v207 we have the following plot. A new paper for the library and also this subcommand is in preparation.

![Alt!](https://github.com/jean-pierreBoth/gsearch/blob/master/GSearch-annembed-GTDBv207.jpg?raw=true)


## Simple case for install

**Pre-built binaries** will be available on release page [binaries](https://github.com/jean-pierreBoth/gsearch/releases/tag/v0.1.3-beta) for major platforms (no need to install but just download and make it executable). We recommend you use the linux one (GSearch_Linux_x86-64_v0.1.3.zip) for linux system in this release page for convenience (only system libraries are required). For macOS, we recommend the binary mac-binaries (GSearch_Darwin_x86-64_v0.1.3.zip or GSearch_Darwin_aarch64_v0.1.3.zip) for corresponding platform (x86-64 or arm64).

## Or if you have conda installed on linux

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gsearch/README.html)

```bash
    ### we suggest python 3.8, so if you do not have one, you can create a new one
    conda create -n python38 python=3.8
    conda activate python38
    conda config --add channels bioconda
    conda install gsearch -c bioconda
    #### a different binary name was used for bioconda channel, you can change it to gsearch later, see below
    gsearchbin -h
    ### change to gsearch binary name
    cp $(which gsearchbin) $CONDA_PREFIX/bin/gsearch
    gsearch -h
```
## or if you are on MacOS and have homebrew installed
```bash
brew tap jianshu93/gsearch
brew install gsearch
gsearch -h

```

Otherwise it is possible to install/compile by yourself (see install section)

```bash

### get the binary for linux (make sure you have recent Linux installed with GCC, e.g., Ubuntu 18.0.4 or above)

wget https://github.com/jean-pierreBoth/gsearch/releases/download/v0.1.3-beta/GSearch_Linux_x86-64_v0.1.3.zip --no-check-certificate
unzip GSearch_Linux_x86-64_v0.1.3.zip

## get the x86-64 binary for macOS
wget https://github.com/jean-pierreBoth/gsearch/releases/download/v0.1.3-beta/GSearch_Darwin_x86-64_v0.1.3.zip --no-check-certificate
unzip GSearch_Darwin_x86-64_v0.1.3.zip
## get the aarch64/arm64 binary for macOS
wget https://github.com/jean-pierreBoth/gsearch/releases/download/v0.1.3-beta/GSearch_Darwin_aarch64_v0.1.3.zip --no-check-certificate
unzip GSearch_Darwin_aarch64_v0.1.3.zip

## get the binary for Windows, ann subcommand is not available for windows for now
wget https://github.com/jean-pierreBoth/gsearch/releases/download/v0.1.3-beta/GSearch_pc-windows-msvc_x86-64_v0.1.3.zip


## Note that for MacOS, you need sudo previlege to allow external binary being executed
* **make it excutable (changed it accordingly on macOS)**
chmod a+x ./gsearch

### put it under your system/usr bin directory (/usr/local/bin/ as an example) where it can be called
mv ./gsearch /usr/local/bin/
### check install
gsearch -h
### check install MacOS, you may need to change the system setup to allow external binary to run by type the following first and use your admin password
sudo spctl --master-disable


```

## usage

```bash
gsearch -h
 ************** initializing logger *****************

Approximate nearest neighbour search for microbial genomes based on MinHash-like metric

Usage: gsearch [OPTIONS] [COMMAND]

Commands:
  tohnsw   Build HNSW graph database from a collection of database genomes based on MinHash-like metric
  add      Add new genome files to a pre-built HNSW graph database
  request  Request nearest neighbors of query genomes against a pre-built HNSW graph database/index
  ann      Approximate Nearest Neighbor Embedding using UMAP-like algorithm
  help     Print this message or the help of the given subcommand(s)

Options:
      --pio <pio>              Parallel IO processing
      --nbthreads <nbthreads>  nb thread for sketching
  -h, --help                   Print help
  -V, --version                Print version

```

We then give here an example of utilization with prebuilt databases.

```bash
### download neighbours for each genomes (fna, fasta, faa et.al. are supported) in query_dir_nt or aa using pre-built database

wget http://enve-omics.ce.gatech.edu/data/public_gsearch/GTDBv207_v2023.tar.gz
tar xzvf ./GTDBv207_v2023.tar.gz

### get test data, we provide 2 genomes at nt, AA and universal gene level

wget https://github.com/jean-pierreBoth/gsearch/releases/download/v0.0.12/test_data.tar.gz --no-check-certificate
tar xzvf ./test_data.tar.gz

cd ./GTDB/nucl
tar -xzvf k16_s12000_n128_ef1600.prob.tar.gz
### request neighbors for nt genomes (here -n is how many neighbors you want to return for each of your query genome)

gsearch --pio 100 --nbthreads 24 request -b ./k16_s12000_n128_ef1600_canonical -r ../../test_data/query_dir_nt -n 50

### or request neighbors for aa genomes (predicted by Prodigal or FragGeneScanRs)

cd ./GTDB/prot
tar xzvf k7_s12000_n128_ef1600.prob.tar.gz
gsearch --pio 100 --nbthreads 24 request -b ./k7_s12000_n128_ef1600_gsearch -r ../../test_data/query_dir_aa -n 50

### or request neighbors for aa universal gene (extracted by hmmer according to hmm files from gtdb, we also provide one in release page)

cd ./GTDB/universal
tar xzvf k5_n128_s1800_ef1600_universal_prob.tar.gz
gsearch --pio 100 --nbthreads 24 request -b ./k5_n128_s1800_ef1600_universal_prob -r ../../test_data/query_dir_universal_aa -n 50

### Building database. database is huge in size, users are welcome to download gtdb database here: (<https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz>) and here (<https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/genomic_files_reps/gtdb_proteins_aa_reps_r207.tar.gz>)

### build database given genome file directory, fna.gz was expected. L for nt and .faa or .faa.gz for --aa. Limit for k is 32 (15 not work due to compression), for s is 65535 (u16) and for n is 255 (u8)

gsearch --pio 2000 --nbthreads 24 tohnsw -d db_dir_nt -s 12000 -k 16 --ef 1600 -n 128 --algo prob
gsearch --pio 2000 --nbthreads 24 tohnsw -d db_dir_aa -s 12000 -k 7 --ef 1600 -n 128 --aa --algo prob

### When there are new genomes  after comparing with the current database (GTDB v207, e.g. ANI < 95% with any genome after searcing, corresponding to >0.875 ProbMinHash distance), those genomes can be added to the database

### old .graph,.data and all .json files will be updated to the new one. Then the new one can be used for requesting as an updated database

gsearch --pio 100 --nbthreads 24 add -b ./k16_s12000_n128_ef1600_canonical -n db_dir_nt (new genomes directory) 

### or add at the amino acid level, in the parameters.json file you can check whether it is DNA or AA data via the "data_t" field
cd ./GTDB/prot
gsearch --pio 100 --nbthreads 24 add -b ./k7_s12000_n128_ef1600_gsearch -n db_dir_nt (new genomes directory in AA format predicted by prodigal/FragGeneScanRs)

### visuzlizing from the tohnsw step at amino acid level (AAI distance), output order of genome files are the same with with seqdict.json
cd ./GTDB/prot
gsearch ann -b ./k7_s12000_n128_ef1600_gsearch --stats --embed
```

## Output explanation

**gsearch.answers.txt** is the default output file in your current directory.  
 For each of your genome in the query_dir, there will be requested N nearest genomes found and sorted by distance (smallest to largest).  
  If one genome in the query does not exist in the output file, meaning at this level (nt or aa), there is no such nearest genomes in the database (or distant away from the best hit in the database), you may then go to amino acid level or universal gene level.

## Dependencies, features and Installation

### Features

- hnsw_rs relies on the crate simdeez to accelerate distance computation. On intel you can build hnsw_rs with the feature simdeez_f

- annembed relies on openblas so you must choose between  the features "annembed_openblas-static" , "annembed_openblas-system" or "annembed_intel-mkl". You may need to install gcc, gfortran and make.
This can be done using the **--features** option as explained below, or by modifying the features section in  Cargo.toml. In that case just fill in the default you want.
- kmerutils provides a feature "withzmq". This feature can be used to store compressed qualities on a server and run requests. It is not necessary in this crate.

## Install/compiling by yourself

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

- On MacOS, which requires dynamic library link (you have to install openblas first):
(note that openblas install lib path is different on M1 MACs).  
So you need to run:

```bash
    brew install openblas
    echo 'export LDFLAGS="-L/usr/local/opt/openblas/lib"' >> ~/.bash_profile
    echo 'export CPPFLAGS="-I/usr/local/opt/openblas/include"' >> ~/.bash_profile
    echo 'export PKG_CONFIG_PATH="/usr/local/opt/openblas/lib/pkgconfig"' >> ~/.bash_profile
    source ~/.bash_profile
    cargo install gsearch --features="annembed_openblas-system"
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
cargo build --release --features="annembed_intel-mkl"
###on MacOS, which requires dynamic library link (install openblas first, see above):
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

- Bacteria/archaea genomes are the newest version of GTDB database (v207, 67,503genomes) (<https://gtdb.ecogenomic.org>), which defines a bacterial speces at 95% ANI. Note that GSearch can also run for even higher resolution species database such as 99% ANI.
- Bacteria/archaea genomes from NCBI/RefSeq until Jan. 2023 (~318,000 genomes)(<https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/>), no clustering at a given ANI threshold.
- Virus data base are based on the JGI IMG/VR database newest version (<https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html>), which also define a virus OTU (vOTU) at 95% ANI.  
- Fungi database are based on the entire RefSeq fungal genomes (retrived via the MycoCosm website), we dereplicated and define a fungal speices at 99.5% ANI.
- All four pre-built databases are available here:<http://enve-omics.ce.gatech.edu/data/gsearch>

## To do list
1. B-bit One Permutation MinHash with Optimal/Fast/BiDirectional densification, which are all locality sensitive hashing scheme and are both more space efficient and faster than classic MinHash with big O O(d+k), see Shrivastava 2017, Mai et.al.,2020 and Jia et.al.,2021. However, the densification step is not mergeable (it introduces randomness), which make it difficult to implement a streaming while hashing & densificaiton scheme. Instead, we must perform hashing & densificaiton when the entire genome/text file is loaded into memory, which requires much more memory in a parallel processing manner.
2. UltraLogLog, Ertl 2023, a significant improvement over HyperLogLog for space/memory for cardinality estimation. We can use inclusion and exclusion rule to estimate Jaccard index (similar to Dashing) after obtaining carinality of two sets despite large variance. A simple but efficient new estimator of cardinality estimation, FGRA (Further Generalized Remaining Area) based on Tau-GRA (Wang and Pettie, 2023) can be used. But we can also map UltraLogLog to corresponding HyperLogLog to use maximum likelihood estimator (MLE) for cardinality (slower), or joint maximum likelihood estimator (JMLE) for intersecion cardinality (then Jaccard index), which are slightly more accurate than FGRA and meet the Cramér-Rao lower bound. However, the Hyperloglog-like scheme is not under the locality sensitive hashing scheme. 

## References

1. Jianshu Zhao, Jean Pierre Both, Luis M. Rodriguez-R and Konstantinos T. Konstantinidis, 2022. GSearch: Ultra-Fast and Scalable Microbial Genome Search by combining Kmer Hashing with Hierarchical Navigable Small World Graphs. *bioRxiv* 2022:2022.2010.2021.513218. [biorxiv](https://www.biorxiv.org/content/10.1101/2022.10.21.513218v2).
