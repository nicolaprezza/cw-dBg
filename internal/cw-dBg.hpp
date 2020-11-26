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

using namespace sdsl;

namespace dBg{

template	<	class dBg_type = int,
				class int_vector = int >
class cw_dBg{

public:


private:

	dBg_type dBg; //the underlying de Bruijn graph
	int_vector deltas; //delta-encoded weights with fast random access

	//marks edges on the dBg that belong to the MST
	rrr_vector<> mst;
	rrr_vector<>::rank_1_type mst_rank;

	//marks sampled nodes on the dBg
	rrr_vector<> sampled;
	rrr_vector<>::rank_1_type sampled_rank;

	//sampled weights
	int_vector samples;

};

}

#endif
