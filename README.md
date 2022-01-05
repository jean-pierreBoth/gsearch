# A rust classifier based on HNSW for prokaryotic genomes

archaea stands for: **A Rust Classifier Hierarchical Navigable SW graphs**

This package (**currently in development**) compute probminhash signature of  bacteria and archaea genomes and stores the id of bacteria and probminhash signature in a Hnsw structure for searching of new request genomes

This package is developped in collaboration with Jianshu Zhao (https://github.com/jianshu93)

# Dependencies and Installation

* pre-compiled binaries are available in the release page for major platforms. For linux based system, no dependencies but system level gfortran must be later than gfortran@5, which means gcc 8.3 or above. If you want to compile from source, see below:


*  Clone hnsw-rs and probminhash or get them from crate.io
    - git clone https://github.com/jean-pierreBoth/hnswlib-rs

    - git clone https://github.com/jean-pierreBoth/probminhash

* Clone kmerutils which is not in crate.io:
  
    - git clone https://github.com/jean-pierreBoth/kmerutils

* Clone ARCHAEA, which is not in crate.io:
    - git clone https://github.com/jean-pierreBoth/archaea
* Clone annembed:
    - git clone https://github.com/jean-pierreBoth/annembed

* Three libraries, zeromq, libsodium and openblas (optional for annembed_f feature) are required to successfully compile. 

```
### if you are using Linux (ubuntu for example), install them first
sudo apt-get install libzmq-dev libsodium-dev openblas

### if you are using MacOS
brew install zeromq
brew install libsodium

### with MacOS monterey, BLAS are installed in the Architecture.framework, but you can still installed it and add library path to environmental variables according to homebrew promt message
brew install openblas

cd archaea
cargo build --release --features annembed_f
## or

cargo build --release


### if you are using anaconda/miniconda3 or on a server, you can install them first, but remember to add library configuration path and dynamic library config path to you environmental variables. Openblas must be installed at system level for MacOS system (static link is not prefered by Apple). Ask your system manager to installed it for you.
conda install -c anaconda zeromq libsodium
## we installed miniconda3 to the home directory
LIBZMQ_LIB_DIR=~/miniconda3/lib LIBZMQ_INCLUDE_DIR=~/miniconda3/include cargo build --release --features annembed_f
### or with out annembed feature, where openblas denpendency is not required
LIBZMQ_LIB_DIR=~/miniconda3/lib LIBZMQ_INCLUDE_DIR=~/miniconda3/include cargo build --release

```

* A dependency is provided as a feature. It uses the crate **annembed** that gives some statistics on the hnsw graph constructed (and will provide some visualization of data).
It is activated by running:
    -   *cargo build --release --features annembed_f*

## build on ARM64/aarch64, rust nightly version only
Nightly rust must be used
```
## install rustup first, and the activate it
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
## setup nightly rust
rustup default nightly
## go to the kmerutils and annembed directory you just cloned, and change the line: hnsw_rs =  {version = "0.1.15"} to hnsw_rs = {path = "../hnswlib-rs"} in both Cargo.toml
## same processure with the above regular compiling.
```
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
  
## Classify
 The classify module is used to assign taxonomy information from requested neighbours to query genomes. Average nucleitide identity will be calculated

## compilation single module
```
cargo build --release --bin tohnsw
cargo build --release --bin request

```
