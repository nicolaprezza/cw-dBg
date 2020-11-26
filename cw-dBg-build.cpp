#include <cassert>
#include <iostream>
#include <chrono>
#include "internal/cw-dBg.hpp"

//using namespace dBg;
using namespace std;
using namespace dbg;

int srate = 50; //sample rate
uint64_t nlines = 0;

format_t format = fastq;

void help(){
	cout << "cw-dBg-build: builds the compressed weighted de Bruijn graph." << endl << endl;
	cout << "Usage: cw-dBg-build [options] <input> <k>" << endl;
	cout << "   Options:"<<endl;
	cout << "   -l <nlines>         use only the first nlines sequences from the input file to build the graph. If set to 0, use all lines. Default: 0."<<endl;
	cout << "   -a                  the input file is fasta. If not specified, it is assumed that the input file is fastq."<<endl;
	cout << "   -s <srate>          sample one out of srate weights. Default: 50."<<endl;
	cout << "   <input>             input fasta/fastq file (see option -a). Mandatory."<<endl;
	cout << "   <k>                 order of the de Bruijn graph. Mandatory."<<endl;
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
	int k = atoi(argv[ptr]);

	cout << "Building compressed weighted de Bruijn graph of input file " << input_file << endl;
	cout << "called as: cw-dBg-build " << (format==fasta?"-a ":"") << "-l " << nlines << " -s " << srate << " " << input_file << " " << k << endl;

	auto t1 = std::chrono::high_resolution_clock::now();

	cw_dBg<> cwdbg(input_file, format);

	auto t2 = std::chrono::high_resolution_clock::now();

	uint64_t elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	cout << "Done. Build time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl;

}
