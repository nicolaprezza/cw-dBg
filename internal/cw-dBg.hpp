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
		return cost_of_weight(dbg.weight_of_edge_(lhs)) > cost_of_weight(dbg.weight_of_edge_(rhs));

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




			// TODO




			if(verbose)
				cout << "Done. Extracting BWT, in/out-degrees, and weights ..." << endl;

			//previous kmer read from kmers.
			__uint128_t prev_kmer = kmers[0];
			out_labels_.push_back(toCHAR(kmers[0] & __uint128_t(7)));
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

					out_labels_.push_back(toCHAR(kmers[i] & __uint128_t(7)));//append to BWT first outgoing edge of this new kmer
					OUT_.push_back(true);//mark that this BWT char is the first of the new kmer

				}else{//same kmer (and thus also (k-1)-mer)

					count++;

					uint8_t curr_char = kmers[i] & __uint128_t(7);
					uint8_t prev_char = prev_kmer & __uint128_t(7);

					//if char of outgoing edge has changed
					if( curr_char != prev_char ){

						//cout << kmer_to_str(kmers[i],k);

						out_labels_.push_back(toCHAR(curr_char));
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
			assert(OUT_.size() == out_labels_.length());
			assert(last.size() == out_labels_.length());

		}//end scope of vector<__uint128_t> kmers;

		MEAN_WEIGHT /= (weights_.size()-padded_kmers);

		C = vector<uint64_t>(5);

		for(auto c : out_labels_) C[toINT(c)]++;

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


		{



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

		auto w = MST_();

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

		return nr_of_nodes - padded_kmers;

	}

	/*
	 * number of edges in the de Bruijn graph.
	 */
	uint64_t number_of_edges(){

		return BWT.size();

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
	 * input: node (represented as an integer) and character
	 * returns: node reached
	 *
	 * if node has no out-edge labeled c, the function returns number_of_nodes()
	 *
	 * c cannot be $
	 *
	 */
	uint64_t move_forward_by_char_(uint64_t node, char c){

		assert(c != '$');

		assert(node <= nr_of_nodes);

		if(node == nr_of_nodes) return node;

		uint64_t idx = start_positions_out_[node];
		uint8_t k = 0;

		while(out_labels_[idx+k] != c and k < out_degree_(node)) ++k;

		return k == out_degree_(node) ? nr_of_nodes : OUT_[idx+k];

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and out_degree_(node)-1
	 * returns: node reached by following the k-th outgoing edge of the node
	 *
	 */
	uint64_t move_forward_by_rank_(uint64_t node, uint8_t k){

		if(node == nr_of_nodes) return node;

		assert(k<out_degree_(node));

		return OUT_[start_positions_out_[node]+k];

	}

	/*
	 * input: edge (represented as pair node, rank of outgoing edge)
	 * returns: position in OUT_ of the edge
	 *
	 */
	uint64_t edge_pos_in_OUT_(edge_t e){

		assert(e.second < out_degree_(e.first));
		return start_positions_out_[e.first] + e.second;

	}

	/*
	 * input: edge (represented as pair node, rank of outgoing edge)
	 * returns: position in array IN (equivalently, in F column) of the edge
	 *
	 */
	uint64_t edge_pos_in_IN_(edge_t e){

		assert(e.second < in_degree_(e.first));

		return start_positions_in_[e.first] + e.second;


	}

	/*
	 * input: node (represented as an integer) and rank between 0 and out_degree_(node)-1
	 * returns: label of k-th outgoing edge of the node
	 *
	 */
	char out_label_(uint64_t node, uint8_t k){

		assert(k<out_degree_(node));

		return out_labels_[start_positions_out_[e.first] + e.second];

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and in_degree_(node)-1
	 * returns: node reached by following the k-th incoming edge of the node
	 *
	 * node cannot be the root 0 ($$$)
	 */
	uint64_t move_backward_(uint64_t node, uint8_t k){

		assert(k<in_degree_(node));
		assert(in_degree_(node)>0);

		assert(node > 0);
		assert(node < number_of_nodes());

		return IN_[start_positions_in_[node]+k];

	}

	/*
	 * returns the in-degree of the node
	 * since this is a dBg and there is only one source node ($$$),
	 * this number is always between 1 and 4
	 *
	 */
	uint8_t in_degree_(uint64_t node){

		assert(node<number_of_nodes());

		return node == nr_of_nodes-1 ? nr_of_nodes - start_positions_in_[node] : start_positions_in_[node+1]-start_positions_in_[node];

	}

	/*
	 * returns the out-degree of the node
	 * this number is always between 0 and 5 (count also $ if present)
	 *
	 */
	uint8_t out_degree_(uint64_t node){

		assert(node<number_of_nodes());

		return node == nr_of_nodes-1 ? nr_of_nodes - start_positions_out_[node] : start_positions_out_[node+1]-start_positions_out_[node];

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

		//mark elements in OUT_ and IN_
		vector<bool> remove_out_edge(OUT_.size(),false);
		vector<bool> remove_in_edge(IN_.size(),false);

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

				for(uint8_t i = 0; i < out_degree_(n) && d<k-1 && out_label_(n,i) != '$' ;++i){

					S.push({move_forward_by_rank_(n,i),d+1});

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
			if(not padded[n] && in_degree_(n)==1 && padded[node = move_backward_(n,0)]){

				while(node != 0){

					//cout << "necessary padded node: " << node << endl;

					assert(padded[node]);
					necessary_node[node] = true;
					assert(in_degree_(node) == 1); //all padded nodes have in-degree 1
					node = move_backward_(node,0);

				}

			}

		}

		assert(necessary_node[0]);

		//update padded_kmers
		padded_kmers = 0;
		for(uint64_t i=0;i<number_of_nodes();++i)
			padded_kmers += padded[i] and necessary_node[i];

		{

			string new_out_labels;
			vector<uint64_t> new_IN;
			vector<uint64_t> new_OUT;
			vector<uint64_t> new_start_positions_in; //for each node, its starting position in IN_
			vector<uint64_t> new_start_positions_out; //for each node, its starting position in OUT_ and out_labels_

			uint64_t start_IN = 0;
			uint64_t start_OUT = 0;

			for(uint64_t n = 0;n<number_of_nodes();++n){

				//append to new_out_labels, new_IN, new_OUT only if the node is necessary
				if(necessary_node[n]){

					new_start_positions_in.push_back(start_IN);
					new_start_positions_out.push_back(start_OUT);

					//in- and out-degree of the node
					auto in_deg = in_degree_(n);
					auto out_deg = out_degree_(n);

					//compute new in-degree
					uint8_t new_in_deg = 0;
					for(uint8_t k = 0;k<in_deg;++k){

						if(n == 0 || necessary_node[move_backward_(n,k)]){

							new_in_deg++;
							new_IN.push_back(n==0?nr_of_nodes:move_backward_(n,k));

						}

					}

					start_IN += new_in_deg;

					assert(new_in_deg>0);

					//compute new out-degree
					uint8_t new_out_deg = 0;

					for(uint8_t k = 0;k<out_deg;++k){

						char c = out_label_(n,k);

						if(c == '$'){

							//insert $ in new_out_labels only if the node has out-degree = 1 (i.e. has only $ as outgoing edge)
							if(out_deg==1){
								new_out_deg++;
								new_OUT.push_back(nr_of_nodes);
								new_out_labels.push_back(c);
							}

						}else{

							//cout << "edge " << char(c) << " node " << n << " out node: " << int(necessary_node[move_forward_by_rank_(n,k)]) << endl;

							if(necessary_node[move_forward_by_rank_(n,k)]){
								new_out_labels.push_back(c);
								new_OUT.push_back(move_forward_by_rank_(n,k));
								new_out_deg++;
							}

						}

					}

					start_OUT += new_out_deg;

					//the root is allowed to have out degree 1 because we mark it as necessary even if it might not be
					//we fix this in the next line
					assert(n == 0 or new_out_deg>0);

					//if root is not actually necessary, we add to it just 1 outgoing edge labeled $
					if(n==0 and new_out_deg == 0){
						new_out_labels.push_back('$');
						new_out_deg=1;
						new_OUT.push_back(nr_of_nodes);
					}

				}

			}

			out_labels_ = new_out_labels;
			IN_ = new_IN;
			OUT_ = new_OUT;
			start_positions_in_ = new_start_positions_in; //for each node, its starting position in IN_
			start_positions_out_ = new_start_positions_out; //for each node, its starting position in OUT_ and out_labels_

			assert(start_positions_in_.size() == start_positions_out_.size());

			nr_of_nodes = start_positions_in_.size();

		}

		//overwirte F_

		F_ = string();
		F_.push_back('$');

		//count number of occurrences of ACGT
		vector<uint64_t> counts(4,0);

		for(auto c:new_out_labels)
			if(c!='$') counts[toINT(c)]++;

		for(uint64_t i=0;i<counts[toINT('A')];++i) F_.push_back('A');
		for(uint64_t i=0;i<counts[toINT('C')];++i) F_.push_back('C');
		for(uint64_t i=0;i<counts[toINT('G')];++i) F_.push_back('G');
		for(uint64_t i=0;i<counts[toINT('T')];++i) F_.push_back('T');

		assert(F_.size() == IN_.size());
		assert(OUT_.size() == out_labels_.size());
		assert();

		//prune weights

		uint64_t idx=0;

		for(uint64_t i=0;i<weights_.size();++i)
			if(necessary_node[i]) weights_[idx++] = weights_[i];

		weights_.resize(number_of_nodes());

		assert(weights_.size() == IN_rank(IN.size()));

	}

	/*
	 * edge is represented as pair <node, rank> where rank is the rank among outgoing edges of node
	 */
	int weight_of_edge_(edge_t e){

		assert(e.first<number_of_nodes());
		assert(e.second < out_degree_(e.first));

		assert(out_label_(e.first, e.second) != '$');

		auto target = OUT_[start_positions_out_[e.first] + e.second];

		assert(target<number_of_nodes());

		return abundance_(e.first) - abundance_(target);

	}

private:

	/*
	 * abundance of node n (represented as integer, i.e. its rank among all nodes)
	 */
	int abundance_(uint64_t n){

		assert(n<number_of_nodes());
		return int(weights_[n]);

	}

	/*
	 * returns character in position i of column F
	 */
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
	 * computes MST forest.
	 * returns weight of the MST
	 */
	uint64_t MST_(){

		parent_in_mst_ = vector<uint64_t>(nr_of_nodes,nr_of_nodes);//roots i will have parent_in_mst_[i] = nr_of_nodes
		mst_out_edges_ = vector<bool>(OUT_.size(),false);

		uint64_t weight = 0; //weight of the MST

		set<uint64_t> not_in_mst;//nodes not yet in the MST

		int n_MST = 0; //number of trees in the MST forest

		//insert all nodes but the root in MST
		for(uint64_t u = 0; u<number_of_nodes(); ++u) not_in_mst.insert(u);

		while(not_in_mst.size()>0){

			n_MST++;

			uint64_t u = *not_in_mst.begin();

			not_in_mst.erase(not_in_mst.find(u));

			//pq contains edges (u,k). Let v = move_forward_by_rank_(u,k). Then edge (u,v) is not in the MST
			priority_queue<edge_t, vector<edge_t>, comp_edge<decltype(*this)> > pq(*this);

			for(uint8_t k = 0; k<out_degree_(u);++k)
				if(out_label_(u,k)!='$')
					pq.push({u,k});

			while(not pq.empty()){

				edge_t e;
				uint64_t v = 0;

				//find a minimum weight edge on the frontier
				do{

					e = pq.top();
					pq.pop();

					v = move_forward_by_rank_(e.first,e.second);

				}while((not pq.empty()) && not_in_mst.find(v) == not_in_mst.end());

				//frontier edge found
				if(not_in_mst.find(v) != not_in_mst.end()){

					weight += cost_of_weight(weight_of_edge_(e));//weight of MST

					not_in_mst.erase(not_in_mst.find(v));
					for(uint8_t k = 0; k<out_degree_(v);++k)
						if(out_label_(v,k)!='$')
							pq.push({v,k});

					parent_in_mst_[v] = e.first;
					mst_out_edges_[start_positions_out_[e.first] + e.second] = true;

				}

			}

		}

		cout << "Number of connected components : " << n_MST << endl;

		assert(not_in_mst.size()==0);

		return weight;

	}

	void build_deltas(){

		//for each node
		for(uint64_t i = 0; i<nr_of_nodes; ++i){

			//if the node is not the root of a tree in the MST forest
			if(parent_in_mst_[i] != nr_of_nodes){

				int w1 = abundance_(parent_in_mst_[i]);
				int w2 = abundance_(i);

				auto encoded_diff = int_to_positive(w1-w2);

				deltas_.push_back(encoded_diff);

			}

		}

	}

	/*
	 * is n a leaf of the MST?
	 */
	bool is_leaf_(uint64_t n){

		uint8_t out_deg = out_degree_(n);

		assert(out_deg>0);

		bool out_edge_found = false;

		uint64_t start = start_positions_out_[n];

		for(uint8_t k=0; k < out_deg;++k){

			out_edge_found = out_edge_found or mst_out_edges_[start+k];

		}

		return not out_edge_found;

	}

	/*
	 * true iff n is a root in the mst forest
	 */
	bool is_root_in_mst_(uint64_t n){

		assert(n < nr_of_nodes);

		return parent_in_mst_[n] != nr_of_nodes;

	}

	/*
	 * marks sampled nodes on the dBg using a tree decomposition that
	 * decomposes the tree in subtrees of size Theta(srate). In the end, marks
	 * the roots of the subtrees using bitvector 'sampled'
	 */
	void macro_tree_decomposition(){

		//rrr_vector<> sampled;
		//rrr_vector<>::rank_1_type sampled_rank;

		sampled_ = vector<bool>(nr_of_nodes,false);

		//post-order visit of the MST
		queue<uint64_t> nodes;

		//size of each component of the tree decomposition
		vector<uint32_t> component_size(nr_of_nodes,1);

		//first insert all leaves
		for(uint64_t i = 0;i<number_of_nodes();++i){

			if(is_leaf_(i))
				nodes.push(i);

		}

		while(not nodes.empty()){

			auto n = nodes.front();
			nodes.pop();

			auto out_deg = out_degree_(n);

			uint16_t sum_of_component_sizes = 1;

			for(uint8_t k=0;k<out_deg;++k){

				uint64_t pos = edge_pos_in_OUT_({n,k});

				//for each MST edge leaving n
				if(mst_out_edges_[pos]){

					sum_of_component_sizes += component_size[move_forward_by_rank_(n,k)];

				}

				if(sum_of_component_sizes < srate){

					component_size[n] = sum_of_component_sizes;

				}else{

					component_size[n] = 0;
					sampled_[n] = 1;

				}

			}

			if(not is_root_in_mst_(n))
				nodes.push(parent_in_mst_[n]);
			else
				sampled_[n] = 1;//roots are always sampled

		}

		assert(sampled_[0] == 1);

		for(uint64_t i=0;i<number_of_nodes();++i){

			if(sampled[i])
				samples_.push_back(weights_[i]);

		}
		cout << samples_.size() << " sampled weights" << endl;

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
	vector<uint32_t> start_positions_out_; //for each node, its starting position in OUT_ and out_labels_

	vector<uint32_t> weights_; //one weight per node

	vector<uint32_t> parent_in_mst_; //one per node. For the roots we write nr_of_nodes
	vector<bool> mst_out_edges_; //marks outgoing edges of the MST, in OUT_ order

	//deltas on MST out edges
	vector<uint64_t> deltas_;

	vector<bool> sampled_; 	//marks sampled nodes on the dBg

	vector<uint64_t> samples_; 	//weight samples

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
