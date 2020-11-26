# cw-dBg: the compressed weighted de Bruijn graph

Author: Nicola Prezza. Joint work with Giuseppe Italiano, Blerina Sinaimeri, Rossano Venturini

### Description

**Warning: experimental code. Runs only on small files (dBg construction is not yet optimized)**

This library builds a compressed representation of the weighted de Bruijn graph. The underlying graph topology is stored using the BOSS representation, while the weights are differentially encoded and sampled on a spanning tree of the graph. Typically, in a 30x-covered dataset the weights are squeezed down to about 3 bits per distinct k-mer. On top of this, the BOSS representation takes about 4 bits per distinct k-mer.

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

to build the compressed weighted de Bruijn graph of order k on the file input (a fastq file by default, or a fasta file if option -a is specified). if option -l nlines is specified, build the graph using only the first nlines sequences from the input file. If option -s srate is specified, sample one out of srate weights (default: srate=50).
