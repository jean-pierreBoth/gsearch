# A Simple database of primitive genoma

This package (**currently in development**) compute probminhash signature of  bacteria and archaea genomes and stores the id of bacteria and probminhash signature in a Hnsw structure

This package is developped in collaboration with Jianshu Zhao

## Sketching of genomes

The sketching and database is done by the module ***tohnsw***.

The sketching of reference genomes can take some time so when done it is stored in 2 structures:
- A Hnsw structure storing rank of data processed and corresponding sketches.
- A Dictionary associating each rank to a fasta id and fasta filename.

The Hnsw structure is dumped *in hnswdump.hnsw.graph* and  *hnswdump.hnsw.data*
The Dictionary is dumped in a json file *seqdict.json*
## Requests

For requests  the module ***request*** is being used. It reloads the dumped files, hnsw and seqdict related
takes a list of fasta files containing requests and for each fasta file dumps the asked number of nearest neighbours.

## compilation

cargo build --release --bin tohnsw