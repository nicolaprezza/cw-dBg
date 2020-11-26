/*
 * dBg.hpp
 *
 *  Created on: Nov 26, 2020
 *      Author: nico
 *
 *
 */

#ifndef INCLUDED_DBG
#define INCLUDED_DBG

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wt_huff.hpp>

using namespace sdsl;
using namespace std;

namespace dbg{

enum format_t {fasta, fastq};

/*
 * bitv_type: used for storing in- and out-degrees
 * str_type: used to store the BWT
 */
template	<	class bitv_type = rrr_vector<>,
				//class bitv_type = bit_vector<>,
				class str_type = wt_huff<>
			>
class dBg{

public:

	dBg(){}

	/*
	 * constructor from fastq/fasta file
	 */
	dBg(string filename, format_t format){

	}

private:

	str_type BWT;

	bitv_type IN;
	typename bitv_type::rank_1_type IN_rank;
	typename bitv_type::select_1_type IN_sel;

	bitv_type OUT;
	typename bitv_type::rank_1_type OUT_rank;
	typename bitv_type::select_1_type OUT_sel;

};

}

#endif
