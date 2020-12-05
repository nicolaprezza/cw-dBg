// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <cassert>
#include <iostream>
#include <chrono>
#include "internal/cw-dBg.hpp"
#include <unordered_map>

//using namespace dBg;
using namespace std;
using namespace dbg;

format_t format = fastq;

uint64_t n = 1000;
int prob = 100;

uint64_t nlines = 0;

void help(){
	cout << "cw-dBg-check: checks a previously constructed weighted de Bruijn graph." << endl << endl;
	cout << "Usage: cw-dBg-check [options] <input_index> <input_fastx>" << endl;
	cout << "   Options:"<<endl;
	cout << "   -l <nlines>         Use only the first nlines sequences from the input file to build the graph. If set to 0, use all lines. Default: 0."<<endl;
	cout << "   -n                  Number of k-mers to test. Default: 1000"<<endl;
	cout << "   -a                  The input file is fasta. If not specified, it is assumed that the input file is fastq."<<endl;
	cout << "   <input_index>       Input index built with cw-dbg-build. Mandatory."<<endl;
	cout << "   <input_fastx>       Fasta/fastq file from which test kmers will be extracted. Must be the same on which the index was built."<<endl;
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

	}else if(s.compare("-l")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -l option." << endl;
			help();
		}

		nlines = atoi(argv[ptr++]);

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


	/*
	 * builds a simple hash table with kmer counts
	 * note: we convert any non-DNA character to A (as it is done in the index)
	 */

	cout << "Building counts with a hash table ... " << endl;

	unordered_map<string, uint64_t> counts;

	{

		ifstream file(input_fastx);

		string str;
		uint64_t n_seq = 0;

		while (std::getline(file, str) and n_seq < nlines) { //getline reads header

			getline(file, str);//getline reads DNA

			for(char & c:str)
				if(c!='A' and c!='C' and c!='G' and c!='T') c = 'A';

			string kmer = str.substr(0,k-1);

			for(int i=k-1;i<str.length();++i){

				kmer += str[i];

				counts[kmer] = counts[kmer]+1;

				kmer = kmer.substr(1,k-1);

			}


			n_seq++;

			if(n_seq%100000==0)
				cout << "processed " << n_seq << " sequences" << endl;

			if(format == fastq){
				getline(file, str);//getline reads +
				getline(file, str);//getline reads quality
			}

		}

	}


	cout << "Testing ... " << endl;

	srand(time(NULL));

	ifstream file(input_fastx);

	uint64_t m = 0;//extracted kmers

	string str;

	int prev_perc = 0;
	int perc = 0;

	while (std::getline(file, str) and m<=n) { //getline reads header

		getline(file, str);//getline reads DNA

		string kmer = str.substr(0,k-1);

		for(int i=k-1;i<str.length() and m<=n;++i){

			kmer += str[i];

			if((rand() % 100) < prob){

				m++;

				if(counts[kmer] != cwdbg[kmer]){

					cout << "Error: " << kmer << " : " << cwdbg[kmer] << " / " << counts[kmer] << endl;
					exit(0);

				}else{

					perc = (m*100)/n;

					if(perc > prev_perc + 4){
						prev_perc = perc;
						cout << (m*100)/n << "% done." << endl;
					}

				}

			}

			kmer = kmer.substr(1,k-1);

		}

		if(format == fastq){
			getline(file, str);//getline reads +
			getline(file, str);//getline reads quality
		}

	}




}
