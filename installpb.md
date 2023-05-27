### Build on MacOS (Intel)

BLAS are installed in the Architecture.framework, but you can still install openblas and add library path to environmental variables according to homebrew promt message. See above annembed part to see how to use a different openblas in intel-mkl library

```bash
### install openblas on intel MACs (note that openblas install lib path is different on M1 MACs)
brew install openblas xz
echo 'export LDFLAGS="-L/usr/local/opt/openblas/lib"' >> ~/.bash_profile
echo 'export CPPFLAGS="-I/usr/local/opt/openblas/include"' >> ~/.bash_profile
echo 'export PKG_CONFIG_PATH="/usr/local/opt/openblas/lib/pkgconfig"' >> ~/.bash_profile
cd gsearch
cargo build --release --features annembed_openblas-system
```
or, if openblas library is not needed

```cargo build --release```

## Build on MacOS ARM64/aarch64

rust nightly version only
Nightly rust must be used
```bash
## install rustup first, and the activate it
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
## setup nightly rust
rustup default nightly
## go to the kmerutils, annembed and gsearch directory you just cloned, and change the line: hnsw_rs =  {version = "0.1.19"} to hnsw_rs = {path = "../hnswlib-rs"} in both Cargo.toml
## same procedure with the above regular compiling.
```
### Crosscompiling for Windows on MacOS or Linux, gcc must be installed (or clang if you are on MacOS via brew)

```
cargo install cargo-xwin
rustup target add x86_64-pc-windows-msvc
cd gsearch
cargo xwin build --target x86_64-pc-windows-msvc

```

### Homology search


The last step involves a homology search using hmmer, which can be directly installed using conda or brew on Intel CPUs. If you are using apple M series ARM64/aarch64 structure, you can have a native support of hmmer folloing the steps:

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
