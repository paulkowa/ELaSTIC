/***
 *  $Id$
 **
 *  File: validate_edges.hpp
 *  Created: May 29, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef VALIDATE_EDGES_HPP
#define VALIDATE_EDGES_HPP

#include <string>
#include <vector>

#include <mpix/data_bucketing.hpp>

#include "SequenceDB.hpp"
#include "config_log.hpp"


inline std::pair<bool, std::string> validate_edges(const AppConfig& opt, AppLog& log, Reporter& report,
						   const std::vector<Sequence>& seqs, SequenceRMA& rma_seq,
						   std::vector<read_pair>& edges,
						   MPI_Comm comm) {
    report << step << "validating edges:" << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    unsigned int n = log.input;

    // generate_edges guarantees that at least one node is local for each edge!

    // divide into remote and local
    unsigned int last =
	(std::partition(edges.begin(), edges.end(), not_local(seqs.front().id, seqs.back().id)) - edges.begin());

    // local edges are easy :-)
    report << info << "processing local edges..." << std::endl;

    sequence_identity ident(seqs, opt.method, opt.kmer, size, rank, n);
    std::transform(edges.begin() + last, edges.end(), edges.begin() + last, ident);

    // here we go with real work
    report << info << "processing remaining edges: " << std::flush;

    id2rank i2r(n, size);
    unsigned int step = (static_cast<double>(last) / 10) + 0.5;

    // sort according to location
    read2rank r2r(rank, i2r);
    std::vector<rank_read_t> rank_id(last);

    if (opt.rma < 2) {
	// order is important!!!
	std::sort(edges.begin(), edges.begin() + last, compare_rank(r2r));
    }

    if (opt.rma == 0) {
	std::transform(edges.begin(), edges.begin() + last, rank_id.begin(), r2r);

	// identify blocks
	std::vector<std::pair<unsigned int, unsigned int> > range;
	unsigned int pos = 0;

	for (unsigned int i = 1; i < last + 1; ++i) {
	    int crank = -1;
	    if (rank_id[pos].first != crank) {
		range.push_back(std::make_pair(pos, i));
		pos = i;
	    }
	} // for i

	// process blocks
	step = (static_cast<double>(range.size()) / 10) + 0.5;
	std::vector<std::string> sv;

	for (unsigned int i = 0; i < range.size(); ++i) {
	    if (i % step == 0) report << "." << std::flush;
	    rma_seq.get(rank_id.begin() + range[i].first, rank_id.begin() + range[i].second, sv);

	    for (unsigned int j = range[i].first; j < range[i].second; ++j) {
		if (i2r(edges[j].id0) == rank) edges[j] = ident(edges[j], sv[j - pos]);
		else edges[j] = ident(sv[j - pos], edges[j]);
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
