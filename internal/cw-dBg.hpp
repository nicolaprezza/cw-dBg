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
#include <internal/dBg.hpp>

using namespace sdsl;
using namespace dbg;
using namespace std;

namespace dbg{

/*
 * cint_vector: used to store the deltas on the MST edges
 */
template	<	class cint_vector = dac_vector_dp<rrr_vector<> >
				//class cint_vector = dac_vector<>
				//class cint_vector = vlc_vector<coder::elias_gamma>
				//class cint_vector = vlc_vector<coder::elias_delta>
			>
class cw_dBg{

public:

	cw_dBg(){}

	/*
	 * constructor from fastq/fasta file
	 */
	cw_dBg(string filename, format_t format){

		G = dBg<>(filename, format);

	}

private:

	dBg<> G; //the underlying de Bruijn graph
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
