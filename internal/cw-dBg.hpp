/*
 * cw-dBg.hpp
 *
 *  Created on: Nov 26, 2020
 *      Author: nico
 *
 *
 */

#ifndef INCLUDED_CWDBG
#define INCLUDED_CWDBG

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <vector>
#include <algorithm>

using namespace sdsl;
using namespace std;

namespace dbg{

enum format_t {fasta, fastq};

/*
 * convert DNA alphabet {$,A,C,G,T} to integers in [0,4] (3 bits per int)
 */
uint8_t toINT(char c){

	switch(c){
		case '$': return 0; break;
		case 'A': case 'a': return 1; break;
		case 'C': case 'c': return 2; break;
		case 'G': case 'g': return 3; break;
		case 'T': case 't': return 4; break;
		default:break;

	}

	return 0;
}

//inverse of the above
char toCHAR(uint8_t x){

	switch(x){
		case 0: return '$'; break;
		case 1: return 'A'; break;
		case 2: return 'C'; break;
		case 3: return 'G'; break;
		case 4: return 'T'; break;
		default:break;

	}

	return '$';

}

/*
 * input: edge (XYZ,W) stored as integer of 128 bits (see "format example" below), character Q stored in 3 bits (function toINT), and order k
 * output: edge (YZW,Q)
 */
__uint128_t edge(__uint128_t kmer, uint8_t c, uint8_t k){

	return (((kmer >> 3) | ((kmer & __uint128_t(7))<<(3*k))) & (~__uint128_t(7))) | c;

}

//has $ symbols in the kmer part?
bool has_dollars(__uint128_t kmer){

	return ((kmer >> 3) &  __uint128_t(7)) == 0;

}

/*
 * bitv_type: used for storing in- and out-degrees (BOSS)
 * str_type: used to store the BWT on alphabet {A,C,G,T,$} (BOSS)
 * cint_vector: used to store the deltas on the MST edges
 */
template	<	class bitv_type = rrr_vector<>,
				//class bitv_type = bit_vector<>,
				class str_type = wt_huff<rrr_vector<> >,
				//class str_type = wt_huff<>
				class cint_vector = dac_vector_dp<rrr_vector<> >
				//class cint_vector = dac_vector<>
				//class cint_vector = vlc_vector<coder::elias_gamma>
				//class cint_vector = vlc_vector<coder::elias_delta>
			>
class cw_dBg{

public:

	cw_dBg(){}

	/*
	 * constructor from fastq/fasta file
	 *
	 * filename: input file path
	 * format: either fasta or fastq
	 * nlines: if 0, process all DNA fragments. Otherwise, only the first nlines
	 * k: de Bruijn graph order (limited to 41)
	 * srate: sample rate
	 */
	cw_dBg(string filename, format_t format, int nlines = 0, uint8_t k = 28, uint16_t srate = 64, bool verbose = true) : k(k), srate(srate){

		assert(k>0 and k<=41);

		string BWT_;
		vector<uint32_t> weights;

		vector<bool> OUT_; //out-degrees: mark first edge (in BWT order) of each new k-mer
		vector<bool> first; //marks first edge (in BWT order) of each new (k-1)-mer

		//begin scope of vector<__uint128_t> kmers;
		{

			ifstream file(filename);

			int read_lines=0;

			/*
			 *  vector storing all (k+1)-mers (i.e. edges) using 3 bits per char
			 *  format example: if k=3 and we have a kmer ACG followed by T, then we store an integer rev(ACG)T = GCAT
			 *
			 */
			vector<__uint128_t> kmers;

			if(verbose)
				cout << "Extracting k-mers from dataset ..." << endl;

			string str;
			while (std::getline(file, str) and (nlines == 0 or read_lines < nlines)) { //getline reads header

				getline(file, str);//getline reads DNA

				// start processing DNA fragment

				//first kmer: ($$$,C), where C is the first letter of str
				__uint128_t kmer = toINT(str[0]);

				kmers.push_back(kmer);

				//push the other kmers
				for(int i=1;i<str.length();++i){

					kmer = edge(kmer, toINT(str[i]),k);
					kmers.push_back(kmer);

				}

				//last kmer: (ACG,$), where ACG was the last kmer in the DNA fragment
				kmer = edge(kmer, toINT('$'),k);
				kmers.push_back(kmer);

				if(format == fastq){
					getline(file, str);//getline reads +
					getline(file, str);//getline reads quality
				}

				read_lines++;

				if(read_lines%100000==0 and verbose)
					cout << "read " << read_lines << " sequences" << endl;

			}

			if(verbose)
				cout << "\nDone. Sorting k-mers ..." << endl;

			sort(kmers.begin(),kmers.end());

			if(verbose)
				cout << "Done. Extracting BWT, in/out-degrees, and weights ..." << endl;

			//previous kmer read from kmers.
			__uint128_t prev_kmer = kmers[0];
			BWT_.push_back(toCHAR(kmers[0] & __uint128_t(7)));
			OUT_.push_back(true);//mark that this BWT char is the first of the new kmer
			first.push_back(true);//mark that this BWT char is the first of the new (k-1)-mer

			uint32_t count = 1;

			for(uint64_t i = 1; i<kmers.size();++i){

				//if (k-1)-mer changes (we will also push in BWT since automatically also the k-mer changes)
				if((kmers[i]>>6) != (prev_kmer>>6)){

					first.push_back(true);

				}

				//if kmer changes: we've already appended edges letter, it's time to append the weight
				if((kmers[i]>>3) != (prev_kmer>>3)){

					//if (k-1)-mer stays the same
					if((kmers[i]>>6) == (prev_kmer>>6)){

						first.push_back(false);

					}

					//we set to 0 the counters of kmers that contain $
					count = has_dollars(prev_kmer)?0:count;
					weights.push_back(count);
					count = 1;

					BWT_.push_back(toCHAR(kmers[i] & __uint128_t(7)));//append to BWT first outgoing edge of this new kmer
					OUT_.push_back(true);//mark that this BWT char is the first of the new kmer

				}else{//same kmer

					count++;

					//if char of outgoing edge has changed
					if( (kmers[i] & __uint128_t(7)) != (prev_kmer & __uint128_t(7)) ){

						BWT_.push_back(toCHAR(kmers[i] & __uint128_t(7)));
						OUT_.push_back(false);//mark that this BWT char is NOT the first of the current kmer
						first.push_back(false);//mark that this BWT char is NOT the first of the current (k-1)-mer

					}

				}

				prev_kmer = kmers[i];

			}

			assert(OUT_.size() == BWT_.length());
			assert(first.size() == BWT_.length());

		}//end scope of vector<__uint128_t> kmers;

		if(verbose)
			cout << "Done. Indexing all structures ..." << endl;





	}

private:

	//parameters:

	uint8_t k; //order
	uint16_t srate; //sample rate

	//BOSS:

	str_type BWT;

	vector<uint64_t> F; //F column

	bitv_type IN;
	typename bitv_type::rank_1_type IN_rank;
	typename bitv_type::select_1_type IN_sel;

	bitv_type OUT;
	typename bitv_type::rank_1_type OUT_rank;
	typename bitv_type::select_1_type OUT_sel;

	//weights:

	cint_vector deltas; //compressed deltas on the edges of the MST

	//marks edges on the dBg that belong to the MST
	rrr_vector<> mst;
	rrr_vector<>::rank_1_type mst_rank;

	//marks sampled nodes on the dBg
	rrr_vector<> sampled;
	rrr_vector<>::rank_1_type sampled_rank;

	//sampled weights
	int_vector<0> samples;

};

}

#endif
