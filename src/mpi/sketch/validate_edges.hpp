/***
 *  $Id$
 **
 *  File: validate_edges.hpp
 *  Created: May 29, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef VALIDATE_EDGES_HPP
#define VALIDATE_EDGES_HPP

#include <string>
#include <vector>

#include <mpix/simple_partition.hpp>

#include "SequenceDB.hpp"
#include "config_log.hpp"


class general_compare {
public:
    explicit general_compare(const AppConfig& opt) : method_(opt.method) {
	switch (opt.method) {
	  case 0:
	      A_ = alignment_identity(opt.gaps[0], opt.gaps[1], opt.gaps[2], opt.gaps[3]);
	      break;
	  case 1:
	      K_ = kmer_identity(opt.kmer, opt.is_dna);
	      break;
	}
    } // general_compare

    unsigned short int operator()(const std::string& s0, const std::string& s1) {
	switch (method_) {
	  case 0: return A_(s0, s1);
	  case 1: return K_(s0, s1);
	}
	return 0;
    } // operator()

private:
    alignment_identity A_;
    kmer_identity K_;

    unsigned int method_;

}; // struct general_compare


inline bool block_compare(const std::pair<unsigned int, unsigned int>& p1, const std::pair<unsigned int, unsigned int>& p2) {
    return (p2.second - p2.first) < (p1.second - p1.first);
} // block_compare


// generate_edges guarantees that at least one node is local for each edge!
inline std::pair<bool, std::string> validate_edges(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
						   const SequenceList& SL, SequenceRMA& rma_seq,
						   std::vector<read_pair>& edges) {
    report << step << "validating edges:" << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // divide into local and remote
    unsigned int mid =
	std::partition(edges.begin(), edges.end(), local(SL.seqs.front().id, SL.seqs.back().id)) - edges.begin();

    // local edges are easy :-)
    report << info << "processing local edges..." << std::endl;

    sequence_identity<general_compare> ident(SL, SL.seqs.front().id, general_compare(opt));
    std::transform(edges.begin(), edges.begin() + mid, edges.begin(), ident);

    // here we go with real work
    report << info << "processing remaining edges: " << std::flush;

    unsigned int n = edges.size() - mid;
    id2rank i2r(SL.N, size);

    if (n > 0) {
	// sort according to location
	read2rank r2r(rank, i2r);
	std::vector<rank_read_t> rank_id(n);

	// sort by rank and target id
	std::sort(edges.begin() + mid, edges.end(), compare_rank(r2r));

	// we transform sorted edges to reads location
	std::transform(edges.begin() + mid, edges.end(), rank_id.begin(), r2r);

	// identify blocks
	const unsigned int SBLOCK = 4096;

	std::vector<std::pair<unsigned int, unsigned int> > edge_range;
	unsigned int pos = 0;

	for (unsigned int i = 1; i < n + 1; ++i) {
	    int crank = -1;
	    if (i < n) crank = rank_id[i].first;
	    if ((rank_id[pos].first != crank) ||
		((rank_id[pos].first == crank) && (rank_id[i].second - rank_id[pos].second > SBLOCK))) {
		edge_range.push_back(std::make_pair(pos, i));
		pos = i;
	    }
	} // for i

	unsigned int step = static_cast<unsigned int>((static_cast<double>(edge_range.size()) / 10) + 0.5);
	if (step == 0) step = 1;

	std::vector<std::string> sv;

	for (unsigned int i = 0; i < edge_range.size(); ++i) {
	    if (i % step == 0) report << "." << std::flush;

	    rma_seq.get(rank_id.begin() + edge_range[i].first, rank_id.begin() + edge_range[i].second, sv);
	    unsigned int start_id = rank_id[edge_range[i].first].second;

	    for (unsigned int j = edge_range[i].first; j < edge_range[i].second; ++j) {
		unsigned int epos = mid + j;
		const std::string& s = sv[rank_id[j].second - start_id];

		if (i2r(edges[epos].id0) == rank) edges[epos] = ident(edges[epos], s);
		else edges[epos] = ident(s, edges[epos]);
	    } // for j
	} // for i

	report << info << std::endl;
    } // if n

    // finalize
    report << step << "cleaning edges..." << std::endl;
    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar(opt.level)), edges.end());

    return std::make_pair(true, "");
} // validate_edges

#endif // VALIDATE_EDGES_HPP
