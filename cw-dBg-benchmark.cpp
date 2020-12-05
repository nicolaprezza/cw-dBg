#include <cassert>
#include <iostream>
#include <chrono>
#include "internal/cw-dBg.hpp"

//using namespace dBg;
using namespace std;
using namespace dbg;

format_t format = fastq;

uint64_t n = 1000;
int prob = 1;


void help(){
	cout << "cw-dBg-benchmark: benchmarks a previously constructed weighted de Bruijn graph." << endl << endl;
	cout << "Usage: cw-dBg-benchmark [options] <input_index> <input_fastx>" << endl;
	cout << "   Options:"<<endl;
	cout << "   -n                  Number of k-mers to test. Default: 1000"<<endl;
	cout << "   -a                  The input file is fasta. If not specified, it is assumed that the input file is fastq."<<endl;
	cout << "   <input_index>       Input index built with cw-dbg-build. Mandatory."<<endl;
	cout << "   <input_fastx>       Fasta/fastq file from which test kmers will be extracted."<<endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-n")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -n option." << endl;
			help();
		}

		n = atoi(argv[ptr++]);

	}else if(s.compare("-a")==0){

		format = fasta;

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

	auto input_index = string(argv[ptr++]);
	auto input_fastx = string(argv[ptr++]);

	cout << "Loading compressed weighted de Bruijn graph from file " << input_index << endl;

	auto t1 = std::chrono::high_resolution_clock::now();

	cw_dBg<> cwdbg(input_index);

	auto t2 = std::chrono::high_resolution_clock::now();

	uint64_t elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	cout << "\nDone. Load time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl;

	cout << "\nDone. Build time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl;

	cout << "Number of kmers " << cwdbg.number_of_distinct_kmers() << endl;
	cout << "Number of dummy nodes " << cwdbg.number_of_padded_kmers() << endl;
	cout << "Total number of nodes (kmers + dummy nodes) " << cwdbg.number_of_nodes() << endl;
	cout << "Number of edges " << cwdbg.number_of_edges() << endl;
	cout << "Max abundance " << cwdbg.max_weight() << endl;
	cout << "Mean abundance " << cwdbg.mean_weight() << endl;

	auto k = cwdbg.get_order();


	srand(time(NULL));

	ifstream file(input_fastx);

	uint64_t m = 0;//extracted kmers

	string str;
	while (std::getline(file, str)) { //getline reads header

		getline(file, str);//getline reads DNA

		string kmer = str.substr(0,k-1);

		for(int i=k-1;i<str.length() and m<=n;++i){

			kmer += str[i];

			if((rand() % 100) < prob){

				cout << kmer << " : " << cwdbg[kmer] << endl;
				m++;

			}

			kmer = kmer.substr(1,k-1);

		}

		if(format == fastq){
			getline(file, str);//getline reads +
			getline(file, str);//getline reads quality
		}

	}




}
