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
#include <stack>
#include <queue>

using namespace sdsl;
using namespace std;

namespace dbg{

typedef pair<uint64_t,uint8_t> edge_t;

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

	return 1;//this includes 'N' characters.
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

uint64_t int_to_positive(int w){

	//encode small integers as small positive integers (negative -> odd, positive -> even)
	return w<0?2*(-w)-1:2*w;

}

//number of bits of x
uint64_t n_bits(uint64_t x){

	 return 64 - __builtin_clzll(x);

}

//gamma-encoding length of the integer w >= 1
uint64_t cost_gamma(uint64_t w){

	assert(w>=1);
	return 2*n_bits(w)-1;

}


//delta-encoding length of the integer w >= 1
uint64_t cost_delta(uint64_t w){

	assert(w>=1);
	return (n_bits(w)-1) + cost_gamma(n_bits(w));

}

//cost function (bits to represent integer). Input: positive integer >= 0
uint64_t cost_of_int(uint64_t w){

	assert(w>=0);
	w += 1; //because gamma and delta can only encode integers >= 1
	return cost_gamma(w);
	//return cost_delta(w);

}


//the cost of a (possibly negative) weight, i.e. the number of bits used to represent it
uint64_t cost_of_weight(int w){

	//return cost_of_int(abs(w));
	return cost_of_int(int_to_positive(w));

}

template<class dbg_type>
class comp_edge{

	public:

	comp_edge(dbg_type& dbg) : dbg(dbg){}

	bool operator() (const edge_t& lhs, const edge_t& rhs) const {

		//minimizes the cost
		//to see the effect of the MST on the overall structure size, turn this > into a <: this way the Maximum Spanning tree will be found
		return cost_of_weight(dbg.weight_of_edge(lhs)) > cost_of_weight(dbg.weight_of_edge(rhs));

	}

	private:

	dbg_type& dbg;

};

/*
 * input: edge (XYZ,W) stored as integer of 128 bits (see "format example" below), character c stored in 3 bits (see function toINT), and order k
 * output: edge (YZW,c)
 */
__uint128_t edge(__uint128_t kmer, uint8_t c, uint8_t k){

	return (((kmer >> 3) | ((kmer & __uint128_t(7))<<(3*k))) & (~__uint128_t(7))) | c;

}

string kmer_to_str(__uint128_t kmer, int k){

	char edge = toCHAR(kmer & __uint128_t(7));

	string km;

	kmer = kmer >> 3;

	for(int i=0;i<k;++i){

		km += toCHAR(kmer & __uint128_t(7));
		kmer = kmer >> 3;

	}

	km += "|";
	km += edge;

	return km;

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
	 * do_not_optimize: turn off space optimization (does not prune dBg)
	 * srate: sample rate
	 */
	cw_dBg(string filename, format_t format, int nlines = 0, uint8_t k = 28, uint16_t srate = 64, bool do_not_optimize = false, bool verbose = true) : k(k), srate(srate){

		assert(k>0 and k<=41);

		string BWT_;
		//vector<uint32_t> weights;

		vector<bool> OUT_; //out-degrees: mark last edge (in BWT order) of each new k-mer
		vector<bool> last; //marks first edge (in BWT order) of each new (k-1)-mer

		if(verbose)
			cout << "Computing how much memory I need to allocate ..." << endl;

		uint64_t pre_allocation = 0;
		uint64_t tot_bases = 0;

		{
			//count how many kmers we will generate (one per base)

			ifstream file(filename);
			int read_lines=0;

			string str;
			while (std::getline(file, str) and (nlines == 0 or read_lines < nlines)) { //getline reads header

				getline(file, str);//getline reads DNA

				pre_allocation += str.length()+1;
				tot_bases += str.length();

				if(format == fastq){
					getline(file, str);//getline reads +
					getline(file, str);//getline reads quality
				}

				read_lines++;

				if(read_lines%1000000==0 and verbose)
					cout << "read " << read_lines << " sequences" << endl;

			}

		}

		//begin scope of vector<__uint128_t> kmers;
		{

			ifstream file(filename);

			if(verbose){
				cout << "Done. Number of bases: " << tot_bases << endl;
				cout << "Trying to allocate " << pre_allocation*16 << " Bytes ... " << endl;
			}

			int read_lines=0;

			/*
			 *  vector storing all (k+1)-mers (i.e. edges) using 3 bits per char
			 *  format example: if k=3 and we have a kmer ACG followed by T, then we store an integer rev(ACG)T = GCAT
			 *
			 */
			vector<__uint128_t> kmers;
			kmers.reserve(pre_allocation);

			if(verbose)
				cout << "Done. Extracting k-mers from dataset ..." << endl;

			string str;
			while (std::getline(file, str) and (nlines == 0 or read_lines < nlines)) { //getline reads header

				getline(file, str);//getline reads DNA

				// start processing DNA fragment

				//first kmer: ($$$,C), where C is the first letter of str

				uint8_t first_char = toINT(str[0]);

				//if(first_char==0) cout << str<<endl;
				assert(first_char!=0);
				__uint128_t kmer = first_char;

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

				if(read_lines%1000000==0 and verbose)
					cout << "read " << read_lines << " sequences" << endl;

			}

			if(verbose)
				cout << "Done. Sorting k-mers ..." << endl;

			sort(kmers.begin(),kmers.end());

			if(verbose)
				cout << "Done. Extracting BWT, in/out-degrees, and weights ..." << endl;

			//previous kmer read from kmers.
			__uint128_t prev_kmer = kmers[0];
			BWT_.push_back(toCHAR(kmers[0] & __uint128_t(7)));
			OUT_.push_back(true);//mark that this BWT char is the first of the new kmer
			last.push_back(true);//mark that this BWT char is the first of the new (k-1)-mer

			uint32_t count = 1;

			int n_nodes = 0;

			//cout << "----" << (n_nodes++) << endl;

			//cout << kmer_to_str(kmers[0],k) << " " << last[0] << endl;

			for(uint64_t i = 1; i<kmers.size();++i){

				//if kmer changes: we've already appended edges letter, it's time to append the weight
				if((kmers[i]>>3) != (prev_kmer>>3)){

					//cout << "----" << (n_nodes++) << endl;

					//cout << kmer_to_str(kmers[i],k);

					//if (k-1)-mer stays the same
					if((kmers[i]>>6) == (prev_kmer>>6)){

						last.push_back(false);

						//cout << " 0" << endl;

					}else{

						last.push_back(true);
						//cout << " 1" << endl;

					}

					//we set to 0 the counters of kmers that contain $
					count = has_dollars(prev_kmer)?0:count;
					weights_.push_back(count); //I'm pushing the weight of the previous kmer

					if(not has_dollars(prev_kmer)){
						MAX_WEIGHT = count>MAX_WEIGHT?count:MAX_WEIGHT;
						MEAN_WEIGHT += count;
					}else{
						padded_kmers++;
					}

					count = 1; //start counting weight of this new kmer

					BWT_.push_back(toCHAR(kmers[i] & __uint128_t(7)));//append to BWT first outgoing edge of this new kmer
					OUT_.push_back(true);//mark that this BWT char is the first of the new kmer

				}else{//same kmer (and thus also (k-1)-mer)

					count++;

					uint8_t curr_char = kmers[i] & __uint128_t(7);
					uint8_t prev_char = prev_kmer & __uint128_t(7);

					//if char of outgoing edge has changed
					if( curr_char != prev_char ){

						//cout << kmer_to_str(kmers[i],k);

						BWT_.push_back(toCHAR(curr_char));
						OUT_.push_back(false);//mark that this BWT char is NOT the first of the current kmer
						last.push_back(false);//mark that this BWT char is NOT the first of the current (k-1)-mer

						//cout << " 0" << endl;

					}

				}

				prev_kmer = kmers[i];

			}

			if(not has_dollars(prev_kmer)){
				MAX_WEIGHT = count>MAX_WEIGHT?count:MAX_WEIGHT;
				MEAN_WEIGHT += count;
			}else{
				padded_kmers++;
			}

			weights_.push_back(count);//push back weight of the lasst kmer
			assert(OUT_.size() == BWT_.length());
			assert(last.size() == BWT_.length());

		}//end scope of vector<__uint128_t> kmers;

		MEAN_WEIGHT /= (weights_.size()-padded_kmers);

		C = vector<uint64_t>(5);

		for(auto c : BWT_) C[toINT(c)]++;

		C[0] = 1;//there is one $ in the F column (it's the only incoming edge of the root kmer $$$)
		C[1] += C[0];
		C[2] += C[1];
		C[3] += C[2];
		C[4] += C[3];

		C[4] = C[3];
		C[3] = C[2];
		C[2] = C[1];
		C[1] = C[0];
		C[0] = 0; //$ starts at position 0 in the F column

		if(verbose)
			cout << "Done. Indexing all structures ..." << endl;

		construct_im(BWT, BWT_.c_str(), 1);

		{

			//flip bits in OUT so that they mark the last (not first) edge of each kmer
			for(uint64_t i=0;i<OUT_.size()-1;++i){

				OUT_[i] = false;
				if(OUT_[i+1]) OUT_[i] = true;

			}
			OUT_[OUT_.size()-1] = true;

			bit_vector out_bv(OUT_.size());
			for(uint64_t i=0;i<OUT_.size();++i)
				out_bv[i] = OUT_[i];

			OUT = bitv_type(out_bv);

		}

		OUT_rank = typename bitv_type::rank_1_type(&OUT);
		OUT_sel = typename bitv_type::select_1_type(&OUT);

		assert(weights_.size() == OUT_rank(OUT.size()));

		{

			//flip bits in last so that they mark the last (not first) edge of each (k-1)-mer
			for(uint64_t i=0;i<last.size()-1;++i){

				last[i] = false;
				if(last[i+1]) last[i] = true;

			}
			last[last.size()-1] = true;


			/*
			 * IN_ will mark with a bit set the last incoming edge of each k-mer
			 * we add 1 because the root kmer $$$ for convention has 1 incoming edge $
			 * (it is the only kmer with incoming edge labeled $)
			 *
			 * to discover how many incoming edges there are, we subtract the number of $
			 * from the BWT size
			 */
			//
			bit_vector IN_(BWT.size() - BWT.rank(BWT.size(),'$')+1, false);

			IN_[0] = true;

			uint64_t cluster_start = 0;//beginning of the current cluster of the (k-1)-mer

			for(uint64_t i=0;i<last.size();++i){

				if(last[i]){//we've found the end of the cluster

					//find previous occurrence of each letter
					uint64_t rn_A = BWT.rank(i+1,'A');
					assert(rn_A<=BWT.rank(BWT.size(),'A'));
					uint64_t last_A = rn_A == 0?BWT.size():BWT.select(rn_A,'A');

					uint64_t rn_C = BWT.rank(i+1,'C');
					assert(rn_C<=BWT.rank(BWT.size(),'C'));
					uint64_t last_C = rn_C == 0?BWT.size():BWT.select(rn_C,'C');

					uint64_t rn_G = BWT.rank(i+1,'G');
					assert(rn_G<=BWT.rank(BWT.size(),'G'));
					uint64_t last_G = rn_G == 0?BWT.size():BWT.select(rn_G,'G');

					uint64_t rn_T = BWT.rank(i+1,'T');
					assert(rn_T<=BWT.rank(BWT.size(),'T'));
					uint64_t last_T = rn_T == 0?BWT.size():BWT.select(rn_T,'T');

					//if the last occurrence of the letter is inside the cluster, then we mark it in the F column
					if(last_A >= cluster_start and last_A <= i){ assert(not IN_[LF(last_A)]); IN_[LF(last_A)] = true;};
					if(last_C >= cluster_start and last_C <= i){ assert(not IN_[LF(last_C)]); IN_[LF(last_C)] = true;};
					if(last_G >= cluster_start and last_G <= i){ assert(not IN_[LF(last_G)]); IN_[LF(last_G)] = true;};
					if(last_T >= cluster_start and last_T <= i){ assert(not IN_[LF(last_T)]); IN_[LF(last_T)] = true;};

					cluster_start = i+1; //in the next position a new cluster starts

				}

			}

			IN = bitv_type(IN_);

		}

		IN_rank = typename bitv_type::rank_1_type(&IN);
		IN_sel = typename bitv_type::select_1_type(&IN);

		//cout << IN_rank(IN.size()) << " " << number_of_nodes() << endl;
		assert(IN_rank(IN.size()) == number_of_nodes());

		//debug only! prints all the data structures
		//print_all();

		if(not do_not_optimize)	prune(verbose);

		//compute MST

		if(verbose)
			cout << "Computing MST ... " << endl;

		auto w = MST();

		if(verbose)
			cout << "Done. Weight of MST: " << w << " bits (" << double(w)/number_of_distinct_kmers() << " bits/kmer if stored with Elias Gamma)" << endl;

		if(verbose)
			cout << "Computing compressed deltas on the MST edges ... " << endl;

		build_deltas();

		if(verbose)
			cout << "Done." << endl;

		if(verbose)
			cout << "Computing MST tree decomposition (samples) ... " << endl;

		macro_tree_decomposition();

		if(verbose)
			cout << "Done." << endl;


		//contruction ended. Free space of weights

		weights_.clear();
		weights_.shrink_to_fit();

	}

	/*
	 * size of the de Bruijn graph
	 */
	uint64_t dbg_size_in_bits(){

		return	8*((size_in_mega_bytes(BWT) +
				size_in_mega_bytes(IN) +
				size_in_mega_bytes(IN_rank) +
				size_in_mega_bytes(IN_sel) +
				size_in_mega_bytes(OUT) +
				size_in_mega_bytes(OUT_rank) +
				size_in_mega_bytes(OUT_sel))*(uint64_t(1)<<20) +
				C.capacity() * 8);

	}

	/*
	 * size (Bytes) of the compressed weights
	 */
	uint64_t weights_size_in_bits(){

		return	8*((size_in_mega_bytes(deltas) +
				size_in_mega_bytes(mst) +
				size_in_mega_bytes(mst_rank) +
				size_in_mega_bytes(sampled) +
				size_in_mega_bytes(sampled_rank) +
				size_in_mega_bytes(samples))*(uint64_t(1)<<20));

	}

	/*
	 * total size (Bytes) of the data structure
	 */
	uint64_t size_in_bits(){

		return dbg_size_in_bits() + weights_size_in_bits();

	}

	/*
	 * number of nodes (distinct kmers + padded kmers) in the de Bruijn graph. Note: also padded kmers are counted.
	 */
	uint64_t number_of_nodes(){

		return nr_of_nodes;

	}

	/*
	 * number of nodes (distinct kmers) in the de Bruijn graph. Note: also padded kmers are counted.
	 */
	uint64_t number_of_padded_kmers(){

		return padded_kmers;
	}

	uint64_t number_of_distinct_kmers(){

		return number_of_nodes() - padded_kmers;

	}

	/*
	 * number of edges in the de Bruijn graph. Dummy edges labeled $ are ignored.
	 */
	uint64_t number_of_edges(){

		return BWT.size() - BWT.rank(BWT.size(),'$');

	}

	uint64_t max_weight(){

		return MAX_WEIGHT;

	}

	double mean_weight(){

		return MEAN_WEIGHT;

	}

	/*
	 * first column of the BWT matrix
	 */
	char F(uint64_t i){

		return toCHAR(F_int(i));

	}

	/*
	 * returns node corresponding to input kmer (as string)
	 *
	 * if kmer does not exist, returns number_of_nodes()
	 *
	 */
	uint64_t find_kmer(string& kmer){

		assert(kmer.length()==k);

		auto range = full_range();

		for(auto c : kmer){

			//cout << "range: " << range.first << " " << range.second << endl;
			range = LF(range,c);

		}

		//cout << "range: " << range.first << " " << range.second << endl;

		if(range.second != range.first+1)
			return number_of_nodes();

		return range.first;

	}

	uint32_t abundance(string& kmer){

		auto idx = find_kmer(kmer);

		return idx == number_of_nodes()?0:weights_[idx]; //TODO: use the compressed weights once 'weights' has been freed

	}

	/*
	 * abundance of node n (represented as integer, i.e. its rank among all nodes)
	 */
	int abundance(uint64_t n){

		assert(n<number_of_nodes());
		return int(weights_[n]); //TODO: use the compressed weights once 'weights' has been freed

	}

	/*
	 * edge is represented as pair <node, rank> where rank is the rank among outgoing edges of node
	 */
	int weight_of_edge(edge_t e){

		assert(e.first<number_of_nodes());
		assert(e.second < out_degree(e.first));

		assert(out_label(e.first, e.second) != '$');

		auto target = move_forward_by_rank(e.first, e.second);

		assert(target<number_of_nodes());

		return abundance(e.first) - abundance(target);

	}

	/*
	 * input: node (represented as an integer) and character
	 * returns: node reached
	 *
	 * if node has no out-edge labeled c, the function returns number_of_nodes()
	 *
	 * c cannot be $
	 *
	 */
	uint64_t move_forward_by_char(uint64_t node, char c){

		assert(c != '$');

		uint64_t nodes = number_of_nodes();

		assert(node <= nodes);

		if(node == nodes) return node;

		//starting point in BWT
		uint64_t start = node == 0 ? 0 : OUT_sel(node)+1;

		while(BWT[start] != c && OUT[start] == 0) start++;

		if(BWT[start] == c)
			node = IN_rank(LF(start));
		else
			node = nodes;

		return node;

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and out_degree(node)-1
	 * returns: node reached by following the k-th outgoing edge of the node
	 *
	 */
	uint64_t move_forward_by_rank(uint64_t node, uint8_t k){

		assert(k<out_degree(node));

		uint64_t pos = node==0?0:OUT_sel(node)+1;

		return IN_rank(LF(pos+k));

	}

	/*
	 * input: edge (represented as pair node, rank of outgoing edge)
	 * returns: position in the BWT of the edge
	 *
	 */
	uint64_t edge_pos_in_bwt(edge_t e){

		assert(e.second < out_degree(e.first));
		return (e.first==0?0:OUT_sel(e.first)+1) + e.second;

	}

	/*
	 * input: edge (represented as pair node, rank of outgoing edge)
	 * returns: position in array IN (equivalently, in F column) of the edge
	 *
	 */
	uint64_t edge_pos_in_IN(edge_t e){

		uint64_t bwt_pos = edge_pos_in_bwt(e);
		return LF(bwt_pos);

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and out_degree(node)-1
	 * returns: label of k-th outgoing edge of the node
	 *
	 */
	char out_label(uint64_t node, uint8_t k){

		assert(k<out_degree(node));

		uint64_t pos = node==0?0:OUT_sel(node)+1;

		return BWT[pos+k];

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and in_degree(node)-1
	 * returns: node reached by following the k-th incoming edge of the node
	 *
	 * node cannot be the root 0 ($$$)
	 */
	uint64_t move_backward(uint64_t node, uint8_t k){

		assert(k<in_degree(node));
		assert(in_degree(node)>0);

		assert(node > 0);
		assert(node < number_of_nodes());
		assert(IN_sel(node)+1 < IN.size());

		uint64_t idx = IN_sel(node)+1+k;

		assert(idx<IN.size());

		return OUT_rank(FL(idx));

	}

	/*
	 * returns the in-degree of the node
	 * since this is a dBg and there is only one source node ($$$),
	 * this number is always between 1 and 4
	 *
	 */
	uint8_t in_degree(uint64_t node){

		assert(node<number_of_nodes());

		return IN_sel(node+1) - (node == 0 ? 0 : IN_sel(node)+1) +1;

	}

	/*
	 * returns the out-degree of the node
	 * this number is always between 0 and 5 (count also $ if present)
	 *
	 */
	uint8_t out_degree(uint64_t node){

		assert(node<number_of_nodes());

		return OUT_sel(node+1) - (node == 0 ? 0 : OUT_sel(node)+1) +1;

	}

	/*
	 * the source of the de Bruijn graph (i.e. node $$$)
	 */
	uint64_t source(){
		return 0;
	}


	/*
	 * removes all the unnecessary padded nodes. The problem with dBgs over read sets is that we add k padded nodes
	 * for each read. Many of these however are not necessary if the dBg is well connected. In particular, only sources
	 * in the dBg (entry points in the connected components) need padded nodes, and those are few.
	 */
	void prune(bool verbose){

		if(verbose){

			cout << "Pruning the dBg ... " << endl;
			cout << "\nStatistics before pruning:" << endl;
			cout << "Number of distinct kmers " << number_of_distinct_kmers() << endl;
			cout << "Number of padded kmers " << number_of_padded_kmers() << endl;
			cout << "Number of nodes (distinct kmers + padded kmers) " << number_of_nodes() << endl;
			cout << "Number of edges " << number_of_edges() << endl;

			cout << "SPACE BEFORE PRUNING: " << endl;
			cout << "  de Bruijn graph (BOSS): " << double(dbg_size_in_bits())/number_of_edges() << " bits per edge" << endl;
			cout << "                          " << double(dbg_size_in_bits())/number_of_nodes() << " bits per node" << endl;
			cout << "                          " << double(dbg_size_in_bits())/number_of_distinct_kmers() << " bits per distinct kmer" << endl;

		}

		//mark all padded nodes
		vector<bool> padded(number_of_nodes(),false);

		//nodes that are necessary (not to be deleted)
		vector<bool> necessary_node(number_of_nodes(),false);

		//mark elements in BWT/OUT and IN
		vector<bool> remove_from_bwt(BWT.size(),false);
		vector<bool> remove_from_in(IN.size(),false);

		{
			//find all nodes that are within distance k-1 from root (i.e. all padded nodes)

			//pairs <node, distance from root>
			stack<pair<uint64_t, uint8_t> > S;
			S.push({0,0}); //push the root

			while(not S.empty()){

				auto p = S.top();
				uint64_t n = p.first;
				uint8_t d = p.second;

				S.pop();

				padded[n] = true;

				for(uint8_t i = 0; i < out_degree(n) && d<k-1 && out_label(n,i) != '$' ;++i){

					S.push({move_forward_by_rank(n,i),d+1});

				}

			}

		}

		//now, for each non-padded kmer X that is preceded only by a padded kmer (i.e. a source of the dBg),
		//mark as necessary all padded kmers that lead from the root to X

		necessary_node[0] = true;

		for(uint64_t n = 0; n<number_of_nodes();++n){

			if(not padded[n]){

				//cout << "necessary non-padded node: " << n << endl;
				necessary_node[n] = true;

			}

			uint64_t node;
			if(not padded[n] && in_degree(n)==1 && padded[node = move_backward(n,0)]){

				while(node != 0){

					//cout << "necessary padded node: " << node << endl;

					assert(padded[node]);
					necessary_node[node] = true;
					assert(in_degree(node) == 1); //all padded nodes have in-degree 1
					node = move_backward(node,0);

				}

			}

		}

		assert(necessary_node[0]);

		//update padded_kmers
		padded_kmers = 0;
		for(uint64_t i=0;i<number_of_nodes();++i)
			padded_kmers += padded[i] and necessary_node[i];

		{

			string newBWT;
			vector<bool> new_IN;
			vector<bool> new_OUT;

			for(uint64_t n = 0;n<number_of_nodes();++n){

				//append to newBWT, newIN, newOUT only if the node is necessary
				if(necessary_node[n]){

					//in- and out-degree of the node
					auto in_deg = in_degree(n);
					auto out_deg = out_degree(n);

					//compute new in-degree
					uint8_t new_in_deg = 0;
					for(uint8_t k = 0;k<in_deg;++k){

						if(n == 0 || necessary_node[move_backward(n,k)])
							new_in_deg++;

					}

					assert(new_in_deg>0);

					for(uint8_t k = 0;k<new_in_deg-1;++k)
						new_IN.push_back(false);

					new_IN.push_back(true);

					//compute new out-degree
					uint8_t new_out_deg = 0;

					for(uint8_t k = 0;k<out_deg;++k){

						char c = out_label(n,k);

						if(c == '$'){

							//insert $ in BWT only if the node has out-degree = 1 (i.e. has only $ as outgoing edge)
							if(out_deg==1){
								new_out_deg++;
								newBWT.push_back(c);
							}

						}else{

							//cout << "edge " << char(c) << " node " << n << " out node: " << int(necessary_node[move_forward_by_rank(n,k)]) << endl;

							if(necessary_node[move_forward_by_rank(n,k)]){
								newBWT.push_back(c);
								new_out_deg++;
							}

						}

					}

					//the root is allowed to have out degree 1 because we mark it as necessary even if it might not be
					//we fix this in the next line
					assert(n == 0 or new_out_deg>0);

					//if root is not actually necessary, we add to it just 1 outgoing edge labeled $
					if(n==0 and new_out_deg == 0){
						newBWT.push_back('$');
						new_out_deg=1;
					}

					for(uint8_t k = 0;k<new_out_deg-1;++k)
						new_OUT.push_back(false);

					new_OUT.push_back(true);

				}

			}

			//re-build BWT
			BWT = str_type();
			construct_im(BWT, newBWT.c_str(), 1);

			bit_vector in_bv(new_IN.size());
			for(uint64_t i=0;i<new_IN.size();++i)
				in_bv[i] = new_IN[i];

			assert(in_bv[in_bv.size()-1]);

			bit_vector out_bv(new_OUT.size());
			for(uint64_t i=0;i<new_OUT.size();++i)
				out_bv[i] = new_OUT[i];

			OUT = bitv_type(out_bv);
			OUT_rank = typename bitv_type::rank_1_type(&OUT);
			OUT_sel = typename bitv_type::select_1_type(&OUT);

			IN = bitv_type(in_bv);
			IN_rank = typename bitv_type::rank_1_type(&IN);
			IN_sel = typename bitv_type::select_1_type(&IN);

		}

		C = vector<uint64_t>(5);

		for(auto c : BWT) C[toINT(c)]++;

		C[0] = 1;//there is one $ in the F column (it's the only incoming edge of the root kmer $$$)
		C[1] += C[0];
		C[2] += C[1];
		C[3] += C[2];
		C[4] += C[3];

		C[4] = C[3];
		C[3] = C[2];
		C[2] = C[1];
		C[1] = C[0];
		C[0] = 0; //$ starts at position 0 in the F column

		//prune weights

		uint64_t idx=0;

		for(uint64_t i=0;i<weights_.size();++i)
			if(necessary_node[i]) weights_[idx++] = weights_[i];

		weights_.resize(number_of_nodes());

		assert(weights_.size() == number_of_nodes());
		assert(OUT.size() == BWT.size());
		assert(IN.size() == BWT.size()-BWT.rank(BWT.size(),'$')+1);
		assert(weights_.size() == IN_rank(IN.size()));

	}

private:

	uint8_t F_int(uint64_t i){

		uint8_t  c = i>=C[1] and i<C[2] ? 1 :
				 i>=C[2] and i<C[3] ? 2 :
				 i>=C[3] and i<C[4] ? 3 :
				 i>=C[4]            ? 4 : 0;

		return c;

	}

	uint64_t LF(uint64_t i){

		char c = BWT[i];

		assert(c!='$');

		return C[toINT(c)] + BWT.rank(i,c);

	}

	pair<uint64_t, uint64_t> full_range(){
		return {0,number_of_nodes()};
	}

	/*
	 * LF function. Input: range OF NODES [begin, end) of nodes reached by path P, and char c
	 * output: nodes reached by path Pc.
	 *
	 * output is an empty range (end <= begin) if input is an empty range or if Pc does not occur
	 *
	 */
	pair<uint64_t, uint64_t> LF(pair<uint64_t, uint64_t> range, char c){

		assert(c!='$');

		uint64_t begin_on_BWT = range.first == 0 ? 0 : OUT_sel(range.first)+1; //inclusive
		uint64_t end_on_BWT = OUT_sel(range.second)+1; //exclusive

		//cout << "begin on BWT: " << begin_on_BWT << " " << end_on_BWT << endl;

		uint64_t c_before_begin = BWT.rank(begin_on_BWT,c);
		uint64_t c_before_end = BWT.rank(end_on_BWT,c);

		//cout << "c before: " << c_before_begin << " " << c_before_end << endl;

		uint64_t first_on_F = C[toINT(c)] + c_before_begin;
		uint64_t last_on_F = C[toINT(c)] + c_before_end;

		//cout << "on F " << first_on_F << " " << last_on_F << endl;

		return {IN_rank(first_on_F), IN_rank(last_on_F)};

	}

	uint64_t FL(uint64_t i){

		assert(i>0);
		assert(i<IN.size());

		uint8_t  c = F_int(i);

		assert(i>=C[c]);

		uint64_t sel = (i-C[c])+1;
		char ch = toCHAR(c);

		assert(sel<=BWT.rank(BWT.size(),ch));
		return BWT.select(sel,ch);

	}

	/*
	 * computes MST forest. Marks in mst edges of IN that are part of the MST
	 * returns weight of the MST
	 */
	uint64_t MST(){

		bit_vector mst_bv(IN.size(),0);

		uint64_t weight = 0; //weight of the MST

		set<uint64_t> not_in_mst;//nodes not yet in the MST

		int n_MST = 0; //number of trees in the MST forest

		//insert all nodes but the root in MST
		for(uint64_t u = 0; u<number_of_nodes(); ++u) not_in_mst.insert(u);

		while(not_in_mst.size()>0){

			n_MST++;

			uint64_t u = *not_in_mst.begin();

			not_in_mst.erase(not_in_mst.find(u));

			//pq contains edges (u,k). Let v = move_forward_by_rank(u,k). Then edge (u,v) is not in the MST
			priority_queue<edge_t, vector<edge_t>, comp_edge<decltype(*this)> > pq(*this);

			for(uint8_t k = 0; k<out_degree(u);++k)
				if(out_label(u,k)!='$')
					pq.push({u,k});

			while(not pq.empty()){

				edge_t e;
				uint64_t v = 0;

				//find a minimum weight edge on the frontier
				do{

					e = pq.top();
					pq.pop();

					v = move_forward_by_rank(e.first,e.second);

				}while((not pq.empty()) && not_in_mst.find(v) == not_in_mst.end());

				//frontier edge found
				if(not_in_mst.find(v) != not_in_mst.end()){

					weight += cost_of_weight(weight_of_edge(e));//weight of MST

					not_in_mst.erase(not_in_mst.find(v));
					for(uint8_t k = 0; k<out_degree(v);++k)
						if(out_label(v,k)!='$')
							pq.push({v,k});

					mst_bv[edge_pos_in_IN(e)] = 1;

				}

			}

		}

		cout << "Number of connected components : " << n_MST << endl;

		assert(not_in_mst.size()==0);

		mst = rrr_vector<>(mst_bv);
		mst_rank = rrr_vector<>::rank_1_type(&mst);

		return weight;

	}

	void build_deltas(){

		vector<uint64_t> deltas_int;

		//for each character in F column
		for(uint64_t i = 0; i<IN.size(); ++i){

			//if the character corresponds to a mst edge
			if(mst[i]){

				//find the two endpoints of the edge
				uint64_t source = OUT_rank(FL(i));
				uint64_t dest = IN_rank(i);

				int w1 = abundance(source);
				int w2 = abundance(dest);

				auto encoded_diff = int_to_positive(w1-w2);

				deltas_int.push_back(encoded_diff);

			}

		}

		assert(mst_rank(mst.size()) == deltas_int.size());

		deltas = cint_vector(deltas_int);

	}

	/*
	 * is n a leaf of the MST?
	 */
	bool is_leaf(uint64_t n){

		uint8_t out_deg = out_degree(n);

		assert(out_deg>0);

		uint64_t first_edge = edge_pos_in_bwt({n,0});
		uint64_t last_edge = first_edge + out_deg-1;

		bool out_edge_found = false;

		while(first_edge<=last_edge){

			if(BWT[first_edge] != '$'){

				auto pos_in_F = LF(first_edge);
				out_edge_found = out_edge_found or mst[pos_in_F];

			}

			first_edge++;

		}

		return not out_edge_found;

	}

	/*
	 * returns parent of n in the MST
	 */
	uint64_t parent_in_mst(uint64_t n){

		assert(not is_root_in_mst(n));

		auto in_deg = in_degree(n);

		uint64_t first_pos = n==0?0:IN_sel(n)+1;
		uint64_t last_pos = first_pos+in_deg-1;

		//there can be only 1 incoming edge in the MST
		assert(mst_rank(last_pos+1) == mst_rank(first_pos)+1);

		while(not mst[first_pos] and first_pos <= last_pos) first_pos++;

		assert(first_pos <= last_pos);

		return OUT_rank(FL(first_pos));

	}

	/*
	 * true iff n is a root in the mst forest
	 */
	bool is_root_in_mst(uint64_t n){

		auto in_deg = in_degree(n);

		uint64_t first_pos = n==0?0:IN_sel(n)+1;
		uint64_t last_pos = first_pos+in_deg-1;

		return mst_rank(last_pos+1) == mst_rank(first_pos);

	}

	/*
	 * true iff position i in the BWT is a MST edge
	 */
	bool mst_edge_on_bwt(uint64_t i){

		if(BWT[i] == '$')
			return false;

		return mst[LF(i)];

	}

	/*
	 * marks sampled nodes on the dBg using a tree decomposition that
	 * decomposes the tree in subtrees of size Theta(srate). In the end, marks
	 * the roots of the subtrees using bitvector 'sampled'
	 */
	void macro_tree_decomposition(){

		//rrr_vector<> sampled;
		//rrr_vector<>::rank_1_type sampled_rank;

		bit_vector sampled_bv(number_of_nodes(),0);

		//post-order visit of the MST
		queue<uint64_t> nodes;

		//size of each component of the tree decomposition
		vector<uint16_t> component_size(number_of_nodes(),1);

		//first insert all leaves
		for(uint64_t i = 0;i<number_of_nodes();++i){

			if(is_leaf(i))
				nodes.push(i);

		}

		while(not nodes.empty()){

			auto n = nodes.front();
			nodes.pop();

			auto out_deg = out_degree(n);

			uint16_t sum_of_component_sizes = 1;

			for(uint8_t k=0;k<out_deg;++k){

				uint64_t pos = edge_pos_in_bwt({n,k});

				//for each MST edge leaving n
				if(mst_edge_on_bwt(pos)){

					sum_of_component_sizes += component_size[move_forward_by_rank(n,k)];

				}

				if(sum_of_component_sizes < srate){

					component_size[n] = sum_of_component_sizes;

				}else{

					component_size[n] = 0;
					sampled_bv[n] = 1;

				}

			}

			if(not is_root_in_mst(n))
				nodes.push(parent_in_mst(n));
			else
				sampled_bv[n] = 1;//roots are always sampled

		}

		assert(sampled_bv[0] == 1);

		sampled = rrr_vector<>(sampled_bv);
		sampled_rank = rrr_vector<>::rank_1_type(&sampled);

		vector<uint32_t> samples_bv;

		for(uint64_t i=0;i<number_of_nodes();++i){

			if(sampled[i])
				samples_bv.push_back(weights_[i]);

		}

		samples = cint_vector(samples_bv);

		cout << samples.size() << " sampled weights" << endl;

	}

	/*
	 * for debug only: prints all structures
	 */
	void print_all(){

		cout << "BOSS: " << endl;

		for(auto b : IN) cout << int(b);
		cout << endl;

		for(int i=0;i<BWT.size();++i)
			cout << F(i);
		cout << endl;

		for(auto c : BWT) cout << c;
		cout << endl;

		for(auto b : OUT) cout << int(b);
		cout << endl;

	}

	//parameters:

	uint8_t k; //order
	uint16_t srate; //sample rate
	uint64_t nr_of_nodes; //number of nodes
	uint64_t MAX_WEIGHT = 0; //max abundance
	double MEAN_WEIGHT = 0; //mean abundance
	uint64_t padded_kmers = 0;//number of nodes corresponding to a k-mers padded with $


	/*
	 * temporary fast data structures used during construction
	 * convention: we end them by underscore
	 */

	//labels of outgoing edges. Corresponds to the BWT, except that there can be unnecessary $ here (will be removed in BWT)
	string out_labels_;
	string F_; //labels of incoming edges
	vector<uint32_t> OUT_; //outgoing edges
	vector<uint32_t> IN_; //incoming edges

	vector<uint32_t> start_positions_in_; //for each node, its starting position in IN_
	vector<uint32_t> start_positions_out_; //for each node, its starting position in OUT_ and BWT_

	vector<uint32_t> weights_; //one weight per node

	vector<uint32_t> parent_in_mst_; //one per node. For the roots we write nr_of_nodes

	/*
	 * compressed data structures
	 */

	//BOSS:

	str_type BWT;

	vector<uint64_t> C; //C array

	bitv_type IN;
	typename bitv_type::rank_1_type IN_rank;
	typename bitv_type::select_1_type IN_sel;

	bitv_type OUT;
	typename bitv_type::rank_1_type OUT_rank;
	typename bitv_type::select_1_type OUT_sel;

	cint_vector deltas; //compressed deltas on the edges of the MST, in IN order (i.e. on the F column)

	//marks edges (in IN order) that belong to the MST
	rrr_vector<> mst;
	rrr_vector<>::rank_1_type mst_rank;

	//marks sampled nodes on the dBg
	rrr_vector<> sampled;
	rrr_vector<>::rank_1_type sampled_rank;

	//sampled weights
	cint_vector samples;

};

}

#endif
