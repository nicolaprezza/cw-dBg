# cw-dBg: the compressed weighted de Bruijn graph

Author: Nicola Prezza. Joint work with Giuseppe Italiano, Blerina Sinaimeri, Rossano Venturini

### Description

**Warning: experimental code. Runs only on small files. dBg construction is not yet optimized and requires 22 Bytes per input base (in RAM) at the moment**

This library builds a compressed representation of the weighted de Bruijn graph. The underlying graph topology is stored using the BOSS representation, while the weights are differentially encoded and sampled on a spanning tree of the graph chosen to minimize the total bit-size of the structure. Results show that on a 20x-covered dataset with 27M of distinct kmers (700Mbases in total) the whole structure takes 5.44 bits per kmer (just 18 MB in total).

### Download

To clone the repository, run:

> git clone http://github.com/nicolaprezza/cw-dBg

### Compile

The library has been tested under linux using gcc 9.2.1. You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

We use cmake to generate the Makefile. Create a build folder in the main cw-dBg folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  cw-dBg-build [-l nlines] [-a] [-s srate] input k

to build the compressed weighted de Bruijn graph of order k on the file input (a fastq file by default, or a fasta file if option -a is specified). if option -l nlines is specified, build the graph using only the first nlines sequences from the input file. If option -s srate is specified, sample one out of srate weights (default: srate=64).
