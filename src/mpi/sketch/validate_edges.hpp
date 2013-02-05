/***
 *  $Id$
 **
 *  File: validate_edges.hpp
 *  Created: May 29, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
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


// generate_edges guarantees that at least one node is local for each edge!
inline std::pair<bool, std::string> validate_edges(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
						   const SequenceList& SL, SequenceRMA& rma_seq,
						   std::vector<read_pair>& edges) {
    report << step << "validating edges:" << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // divide into remote and local
    unsigned int last =
	std::partition(edges.begin(), edges.end(), not_local(SL.seqs.front().id, SL.seqs.back().id)) - edges.begin();

    // local edges are easy :-)
    report << info << "processing local edges..." << std::endl;

    sequence_identity<general_compare> ident(SL, SL.seqs.front().id, general_compare(opt));
    std::transform(edges.begin() + last, edges.end(), edges.begin() + last, ident);

    // here we go with real work
    report << info << "processing remaining edges: " << std::flush;

    id2rank i2r(SL.N, size);

    unsigned int step = (static_cast<double>(last) / 10) + 0.5;
    if (step == 0) step = 1;

    // sort according to location
    read2rank r2r(rank, i2r);
    std::vector<rank_read_t> rank_id(last);

    // rma0 - sequences prefetched in blocks, blocks in random order
    // rma1 - sequences accessed one-by-one, in ordered fashion
    // rma2 - sequences accessed one-by-one, in random order

    // sort by rank and local id
    if (opt.rma < 2) {
	std::sort(edges.begin(), edges.begin() + last, compare_rank(r2r));
    }

    if (opt.rma == 0) {
	// we transform sorted edges to reads location
	std::transform(edges.begin(), edges.begin() + last, rank_id.begin(), r2r);

	// identify blocks
	std::vector<std::pair<unsigned int, unsigned int> > range;
	unsigned int pos = 0;

	for (unsigned int i = 1; i < last + 1; ++i) {
	    int crank = -1;
	    if (i < last) crank = rank_id[i].first;
	    if (rank_id[pos].first != crank) {
		range.push_back(std::make_pair(pos, i));
		pos = i;
	    }
	} // for i

	std::srand(std::time(0));
	std::random_shuffle(range.begin(), range.end());

	// process blocks
	step = (static_cast<double>(range.size()) / 10) + 0.5;
	if (step == 0) step = 1;

	std::vector<std::string> sv;

	for (unsigned int i = 0; i < range.size(); ++i) {
	    if (i % step == 0) report << "." << std::flush;
	    rma_seq.get(rank_id.begin() + range[i].first, rank_id.begin() + range[i].second, sv);

	    for (unsigned int j = range[i].first; j < range[i].second; ++j) {
		if (i2r(edges[j].id0) == rank) edges[j] = ident(edges[j], sv[j - range[i].first]);
		else edges[j] = ident(sv[j - range[i].first], edges[j]);
	    } // for j

	} // for i
    } else {
	std::string s = "";

	for (unsigned int i = 0; i < last; ++i) {
	    if (i % step == 0) report << "." << std::flush;
	    if (i2r(edges[i].id0) == rank) {
		// id0 is local
		s = rma_seq.get(edges[i].id1);
		edges[i] = ident(edges[i], s);
	    } else {
		s = rma_seq.get(edges[i].id0);
		edges[i] = ident(s, edges[i]);
	    }
	} // for i
    } // if opt.rma

    report << "" << std::endl;

    report << "cleaning edges..." << std::endl;
    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar(opt.level)), edges.end());

    return std::make_pair(true, "");
} // validate_edges

#endif // VALIDATE_EDGES_HPP
