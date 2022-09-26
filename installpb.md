## installation help with libzmq

###  zmq install on linux server without sudo privilege (install miniconda3 first):

* Install via conda (recommended):

  `conda activate`

  `conda install zeromq`

  

  Install Rustup.  On linux: `conda install -c milesgranger rustup` or on MacOs :  
  `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

  
  `rustup default nightly`
  
  change to you miniconda installation path
  
  `a=$(which conda)`
  
  `LIBZMQ_LIB_DIR=${a%/*/*}/lib LIBZMQ_INCLUDE_DIR=${a%/*/*}/include cargo install archaea`
  
  `cargo install --git https://gitlab.com/Jianshu_Zhao/fraggenescanrs`
  
  `conda install hmmer`


 

### zmq and libsodium install on linux

```
sudo apt-get install libzmq-dev libsodium-dev openblas
```

### zmq and libsodium on MacOS

```
brew install zeromq  
brew install libsodium
```

### with MacOS monterey

BLAS are installed in the Architecture.framework, but you can still installed it and add library path to environmental variables according to homebrew promt message. See above annembed part to see how to use a different openblas in intel-mkl library
brew install openblas

```
cd archaea
cargo build --release --features annembed_f  (if annembed is needed or cargo build --release=
```
or, if openblas library is not needed

```cargo build --release```


## using anaconda/miniconda3 or on a server

you can install them first, but remember to add library configuration path and dynamic library config path to you environmental variables. Openblas must be installed at system level for MacOS system (static link is not prefered by Apple). Ask your system manager to install it for you.
### with installed miniconda3 to the home directory
```
LIBZMQ_LIB_DIR=~/miniconda3/lib LIBZMQ_INCLUDE_DIR=~/miniconda3/include cargo build --release --features annembed_f
```

or without annembed feature, in this case openblas denpendency is not required:
```
LIBZMQ_LIB_DIR=~/miniconda3/lib LIBZMQ_INCLUDE_DIR=~/miniconda3/include cargo build --release
```




## Build on ARM64/aarch64,

rust nightly version only
Nightly rust must be used
```bash
## install rustup first, and the activate it
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
## setup nightly rust
rustup default nightly
## go to the kmerutils, annembed and archaea directory you just cloned, and change the line: hnsw_rs =  {version = "0.1.15"} to hnsw_rs = {path = "../hnswlib-rs"} in both Cargo.toml
## same procedure with the above regular compiling.
```


