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

#include "../StaticWSQueue.hpp"

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "create_smatrix.hpp"

#include "iomanip.hpp"


class general_compare {
public:
    explicit general_compare(const AppConfig& opt) : method_(opt.method) {
	bio::scoring_matrix sm;
	int g, h;

	if (opt.method > 0) create_smatrix(opt.gaps, opt.is_dna, sm, g, h);

	switch (opt.method) {
	  case 0:
	      Kid_ = kmer_identity(opt.kmer, opt.is_dna);
	      break;

	  case 1:
	      cAcdhit_ = cfe_alignment_cdhit_identity(sm, g, h);
	      break;

	  case 2:
	      Acdhit_ = alignment_cdhit_identity(sm, g, h);
	      break;

	  case 3:
	      cAblast_ = cfe_alignment_blast_identity(sm, g, h);
	      break;

	  case 4:
	      Ablast_ = alignment_blast_identity(sm, g, h);
	      break;
	} // switch
    } // general_compare

    unsigned short int operator()(const std::string& s0, const std::string& s1) {
	const std::string* sa = &s0;
	const std::string* sb = &s1;

	if (s1.size() < s0.size()) std::swap(sa, sb);

	switch (method_) {
	  case 0: return Kid_(*sa, *sb);
	  case 1: return Acdhit_(*sa, *sb);
	  case 2: return cAcdhit_(*sa, *sb);
	  case 3: return Ablast_(*sa, *sb);
	  case 4: return cAblast_(*sa, *sb);
	}

	return 0;
    } // operator()

private:
    kmer_identity Kid_;
    alignment_cdhit_identity Acdhit_;
    cfe_alignment_cdhit_identity cAcdhit_;
    alignment_blast_identity Ablast_;
    cfe_alignment_blast_identity cAblast_;

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

    if (n > 0) {
	// sort according to location
	id2rank i2r(SL.N, size);
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
    report << info << "cleaning edges..." << std::endl;
    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar(opt.level)), edges.end());

    return std::make_pair(true, "");
} // validate_edges


inline std::pair<bool, std::string> validate_edges_ws(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
						      const SequenceList& SL, SequenceRMA& rma_seq,
						      std::vector<read_pair>& edges) {
    report << step << "validating edges:" << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    report << info << "creating tasks..." << std::endl;

    // divide into local and remote
    unsigned int mid =
	std::partition(edges.begin(), edges.end(), local(SL.seqs.front().id, SL.seqs.back().id)) - edges.begin();

    // sort remote according to location
    id2rank i2r(SL.N, size);
    read2rank r2r(rank, i2r);

    std::sort(edges.begin() + mid, edges.end(), compare_rank(r2r));

    // transform remote edges to reads location
    unsigned int n = edges.size() - mid;
    std::vector<rank_read_t> rank_id(n);

    std::transform(edges.begin() + mid, edges.end(), rank_id.begin(), r2r);

    // get tasks list
    typedef StaticWSQueue<read_pair> ws_queue_type;
    const unsigned int SBLOCK = 4096;

    std::vector<ws_queue_type::range_type> tasks;
    unsigned int pos = 0;

    for (unsigned int i = 1; i < n + 1; ++i) {
	int crank = -1;
	if (i < n) crank = rank_id[i].first;
	if ((rank_id[pos].first != crank) ||
	    ((rank_id[pos].first == crank) && (rank_id[i].second - rank_id[pos].second > SBLOCK))) {
	    tasks.push_back(ws_queue_type::make_range(pos, i));
	    pos = i;
	}
    } // for i

    report << info << "initializing queue..." << std::endl;

    // create work queue from remote edges
    ws_queue_type wsq(comm);

    if (wsq.init(edges.begin() + mid, edges.end(), tasks.begin(), tasks.end()) == false) {
	return std::make_pair(false, "task queue failed to initialize");
    }

    // process local edges
    report << info << "processing local edges..." << std::endl;

    sequence_identity<general_compare> ident(SL, SL.seqs.front().id, general_compare(opt));
    edges.resize(mid);

    std::transform(edges.begin(), edges.end(), edges.begin(), ident);
    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar(opt.level)), edges.end());

    // process tasks
    report << info << "processing remaining edges, be patient..." << std::endl;

    const read_pair* first = 0;
    const read_pair* last = 0;

    // here one read is local
    std::vector<std::string> sv;

    while (wsq.get(first, last) == true) {
	unsigned int l = last - first;

	rank_id.resize(l);
	std::transform(first, last, rank_id.begin(), r2r);

	rma_seq.get(rank_id.begin(), rank_id.end(), sv);
	unsigned int start_id = rank_id[0].second;

	for (unsigned int j = 0; j < l; ++j) {
	    const std::string& s = sv[rank_id[j].second - start_id];
	    const read_pair& e = first[j];

	    if (i2r(e.id0) == rank) edges.push_back(ident(e, s));
	    else edges.push_back(ident(s, e));
	} // for j

	wsq.progress();
    } // while wsq.get

    // we start stealing process
    report << info << "almost done..." << std::endl;

    int vrank = -1;
    read_pair* sfirst = 0;
    read_pair* slast = 0;

    std::string s0;

    while (wsq.steal(vrank, sfirst, slast) == true) {
	unsigned int l = slast - sfirst;
	read2rank vr2r(vrank, i2r);

	rank_id.resize(l);
	std::transform(sfirst, slast, rank_id.begin(), vr2r);

	rma_seq.get(rank_id.begin(), rank_id.end(), sv);
	unsigned int start_id = rank_id[0].second;

	for (unsigned int j = 0; j < l; ++j) {
	    const std::string& s = sv[rank_id[j].second - start_id];
	    const read_pair& e = sfirst[j];

	    if (i2r(e.id0) == rank_id[j].first) {
		s0 = rma_seq.get(e.id1);
		edges.push_back(ident(e, s, s0));
	    } else {
		s0 = rma_seq.get(e.id0);
		edges.push_back(ident(e, s0, s));
	    }
	} // for

	delete[] sfirst;
    } // while wsq.steal

    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar(opt.level)), edges.end());
    wsq.finalize();

    return std::make_pair(true, "");
} // validate_edges_ws

#endif // VALIDATE_EDGES_HPP
