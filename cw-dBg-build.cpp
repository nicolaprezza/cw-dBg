#include <cassert>
#include <iostream>
#include <chrono>
#include "internal/cw-dBg.hpp"

//using namespace dBg;
using namespace std;
using namespace dbg;

uint16_t srate = 64; //sample rate
uint64_t nlines = 0;

format_t format = fastq;

bool pause_ = false;
bool do_not_optimize = false;

void help(){
	cout << "cw-dBg-build: builds the compressed weighted de Bruijn graph." << endl << endl;
	cout << "Usage: cw-dBg-build [options] <input> <k>" << endl;
	cout << "   Options:"<<endl;
	cout << "   -l <nlines>         Use only the first nlines sequences from the input file to build the graph. If set to 0, use all lines. Default: 0."<<endl;
	cout << "   -s <srate>          Sample one out of srate weights. Default: 50."<<endl;
	cout << "   -a                  The input file is fasta. If not specified, it is assumed that the input file is fastq."<<endl;
	cout << "   -o                  Turn off space optimization (does not prune the dBg)."<<endl;
	cout << "   -p                  Pause exectution before and after construction in order to allow measuring RAM."<<endl;
	cout << "   <input>             Input fasta/fastq file (see option -a). Mandatory."<<endl;
	cout << "   <k>                 Order of the de Bruijn graph in [1,41]. Mandatory."<<endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-l")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -l option." << endl;
			help();
		}

		nlines = atoi(argv[ptr++]);

	}else if(s.compare("-a")==0){

		format = fasta;

	}else if(s.compare("-p")==0){

		pause_ = true;

	}else if(s.compare("-o")==0){

		do_not_optimize = true;

	}else if(s.compare("-s")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -s option." << endl;
			help();
		}

		srate = atoi(argv[ptr++]);

	}else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

	//parse options

	int ptr = 1;

	if(argc<3) help();

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	auto input_file = string(argv[ptr++]);
	uint8_t k = atoi(argv[ptr]);

	if(k>41 or k==0){
		cout << "Error: k must be in [1,41]" << endl;
		help();
	}

	cout << "Building compressed weighted de Bruijn graph of input file " << input_file << endl;
	cout << "called as: cw-dBg-build " << (format==fasta?"-a ":"") << "-l " << nlines << " -s " << srate << " " << input_file << " " << int(k) << endl;

	if(pause_){

		cout << "Paused: construction has to start. Press enter to continue" << endl;
		cin.ignore();

	}

	auto t1 = std::chrono::high_resolution_clock::now();

	//cw_dBg<bit_vector, wt_huff<> > cwdbg(input_file, format, nlines, k, srate, true); //fast - uses uncompressed vectors
	cw_dBg<> cwdbg(input_file, format, nlines, k, srate, do_not_optimize, true); //slow but very small - uses rrr-compressed bit-vectors everywhere

	string km;
	km = "CGA";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "GAC";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "TAC";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "GTC";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "ACG";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "TCG";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "ACT";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "CGT";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;
	km = "AAA";	cout << km << " " << cwdbg.find_kmer(km) << " " << cwdbg.abundance(km) << endl;

	//for(int i=0;i<11;++i) cout << i << " " << int(cwdbg.in_degree(i)) << " " << int(cwdbg.out_degree(i)) << endl;

	auto t2 = std::chrono::high_resolution_clock::now();

	if(pause_){

		cout << "Paused: construction has ended. Press enter to continue" << endl;
		cin.ignore();

	}

	uint64_t elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	cout << "Done. Build time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl << endl;

	cout << "Number of distinct kmers " << cwdbg.number_of_distinct_kmers() << endl;
	cout << "Number of padded kmers " << cwdbg.number_of_padded_kmers() << endl;
	cout << "Number of nodes (distinct kmers + padded kmers) " << cwdbg.number_of_nodes() << endl;
	cout << "Number of edges " << cwdbg.number_of_edges() << endl;
	cout << "Max abundance " << cwdbg.max_weight() << endl;
	cout << "Mean abundance (only on distinct kmers) " << cwdbg.mean_weight() << endl;

	cout << "SPACE: " << endl;
	cout << "  de Bruijn graph (BOSS): " << double(cwdbg.dbg_size_in_bits())/cwdbg.number_of_edges() << " bits per edge" << endl;
	cout << "                          " << double(cwdbg.dbg_size_in_bits())/cwdbg.number_of_nodes() << " bits per node" << endl;
	cout << "                          " << double(cwdbg.dbg_size_in_bits())/cwdbg.number_of_distinct_kmers() << " bits per distinct kmer" << endl;
	cout << "  compressed weights:     " << double(cwdbg.weights_size_in_bits())/cwdbg.number_of_edges() << " bits per edge" << endl;
	cout << "                          " << double(cwdbg.weights_size_in_bits())/cwdbg.number_of_nodes() << " bits per node" << endl;
	cout << "                          " << double(cwdbg.weights_size_in_bits())/cwdbg.number_of_distinct_kmers() << " bits per distinct kmer" << endl;
	cout << "  total:                  " << double(cwdbg.size_in_bits())/cwdbg.number_of_edges() << " bits per edge" << endl;
	cout << "                          " << double(cwdbg.size_in_bits())/cwdbg.number_of_nodes() << " bits per node" << endl;
	cout << "                          " << double(cwdbg.size_in_bits())/cwdbg.number_of_distinct_kmers() << " bits per distinct kmer" << endl<<endl;




}
