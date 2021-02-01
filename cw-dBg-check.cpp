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

int q = 1000000;

bool check = false;

int prob = 20; //probability of extracting a kmer

uint64_t nlines = 0;

void help(){
	cout << "cw-dBg-check: checks a previously constructed weighted de Bruijn graph." << endl << endl;
	cout << "Usage: cw-dBg-check [options] <input_index> <input_fastx>" << endl;
	cout << "   Options:"<<endl;
	cout << "   -q <arg>            Extract and test the structure on the first maximum <arg> k-mers in the dataset. Default: " << q << "." <<endl;
	cout << "   -a                  The input file is fasta. If not specified, it is assumed that the input file is fastq."<<endl;
	cout << "   -c                  Check correctness of the structure against a classic hash (space-consuming!!). Default: false."<<endl;
	cout << "   <input_index>       Input index built with cw-dbg-build. Mandatory."<<endl;
	cout << "   <input_fastx>       Fasta/fastq file from which test kmers will be extracted. Must be the same on which the index was built."<<endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-a")==0){

		format = fasta;

	}/*else if(s.compare("-l")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -l option." << endl;
			help();
		}

		nlines = atoi(argv[ptr++]);

	}*/
	else if(s.compare("-c")==0){

		check = true;

	}else if(s.compare("-q")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -lq option." << endl;
			help();
		}

		q = atoi(argv[ptr++]);

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

	cout << "k = " << int(cwdbg.get_order()) << endl;
	cout << "Number of kmers " << cwdbg.number_of_distinct_kmers() << endl;
	cout << "Number of dummy nodes " << cwdbg.number_of_padded_kmers() << endl;
	cout << "Total number of nodes (kmers + dummy nodes) " << cwdbg.number_of_nodes() << endl;
	cout << "Number of edges " << cwdbg.number_of_edges() << endl;
	cout << "Max abundance " << cwdbg.max_weight() << endl;
	cout << "Mean abundance " << cwdbg.mean_weight() << endl;

	auto k = cwdbg.get_order();


	vector<string> kmers;
	kmers.reserve(q);


	cout << "Extracting kmers ... " << endl;

	srand(time(NULL));

	ifstream file(input_fastx);

	uint64_t m = 0;//extracted kmers

	string str;

	int prev_perc = 0;
	int perc = 0;

	while (std::getline(file, str) and m<q) { //getline reads header

		getline(file, str);//getline reads DNA

		for(char & c:str)
			if(c!='A' and c!='C' and c!='G' and c!='T') c = 'A';

		string kmer = str.substr(0,k-1);

		for(int i=k-1;i<str.length() and m<q;++i){

			kmer += str[i];

			if((rand() % 100) < prob){

				m++;

				kmers.push_back(kmer);

				////process kmer

			}

			kmer = kmer.substr(1,k-1);

		}

		if(format == fastq){
			getline(file, str);//getline reads +
			getline(file, str);//getline reads quality
		}

	}

	cout << "Done." << endl;


	if(check){

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

			__uint128_t kmer = 0;

			while (std::getline(file, str) and (nlines == 0 or n_seq < nlines)) { //getline reads header

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

		cout << "Done." << endl;

		cout << "Testing ... " << endl;

		prev_perc = 0;
		perc = 0;

		m=0;

		for(auto & kmer : kmers){

			m++;

			if(counts[kmer] != cwdbg[kmer]){

				cout << "Error: " << kmer << " : " << cwdbg[kmer] << " / " << counts[kmer] << endl;
				exit(0);

			}else{

				perc = (m*100)/kmers.size();

				if(perc > prev_perc + 4){
					prev_perc = perc;
					cout << (m*100)/kmers.size() << "% done." << endl;
				}

			}
		}

	}


	cout << "Benchmarking the structure on " << kmers.size() << " abundance queries." << endl;

	prev_perc = 0;
	perc = 0;

	auto tb1 = std::chrono::high_resolution_clock::now();

	m=0;

	for(auto & kmer : kmers){

		if(cwdbg[kmer] != 0){

			m++;

		}

		perc = (m*100)/kmers.size();

		if(perc > prev_perc + 4){
			prev_perc = perc;
			cout << (m*100)/kmers.size() << "% done." << endl;
		}

	}

	auto tb2 = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(tb2 - tb1).count();
	cout << "\nDone. " << double(elapsed)/kmers.size() << " microseconds/query" << endl;


}
