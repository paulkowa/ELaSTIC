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

#include "Sequence.hpp"
#include "config_log.hpp"


inline std::pair<bool, std::string> validate_edges(const AppConfig& opt, AppLog& log, Reporter& report,
						   const std::vector<Sequence>& seqs,
						   std::vector<read_pair>& edges,
						   MPI_Comm comm) {
    report << step << "validating edges..." << std::endl;
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    unsigned int mid = 0;
    unsigned int last = 0;

    last = std::partition(edges.begin(), edges.end(), not_local(seqs.front().id, seqs.back().id)) - edges.begin();
    mid = std::partition(edges.begin(), edges.begin() + last, not_collocated(log.input, size)) - edges.begin();

    unsigned int loc[3];
    loc[0] = edges.size() - last;
    loc[1] = last - mid;
    loc[2] = mid;

    unsigned int min[3] = { 0, 0, 0 };
    unsigned int max[3] = { 0, 0, 0 };
    unsigned int tot[3] = { 0, 0, 0 };

    MPI_Reduce(loc, min, 3, MPI_UNSIGNED, MPI_MIN, 0, comm);
    MPI_Reduce(loc, max, 3, MPI_UNSIGNED, MPI_MAX, 0, comm);

    report << info << "local edges distribution: [" << min[0] << "," << max[0] << "]" << std::endl;
    report << info << "collocated edges distribution: [" << min[1] << "," << max[1] << "]" << std::endl;
    report << info << "remaining edges distribution: [" << min[2] << "," << max[2] << "]" << std::endl;

    // first local edges (i.e. reads and edge are available)
    report << info << "processing local edges..." << std::endl;

    std::transform(edges.begin() + last, edges.end(), edges.begin() + last, identity(seqs, opt.method, opt.kmer));
    edges.erase(std::remove_if(edges.begin() + last, edges.end(), not_similar(opt.level)), edges.end());

    // now we group edges on the processors that can process them (edges collocated)
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    read_pair* efirst = 0;
    read_pair* elast = 0;

    boost::tie(efirst, elast) = mpix::data_bucketing(edges.begin() + mid, edges.begin() + last, read_pair_local_hash(log.input, size),
						     MPI_READ_PAIR, 0, comm);
    edges.erase(edges.begin() + mid, edges.begin() + last);

    // we process received batch
    report << info << "processing collocated edges..." << std::endl;

    std::transform(efirst, elast, efirst, identity(seqs, opt.method, opt.kmer));
    elast = std::remove_if(efirst, elast, not_similar(opt.level));
    std::copy(efirst, elast, std::back_inserter(edges));

    delete[] efirst;

    MPI_Type_free(&MPI_READ_PAIR);

    unsigned int per = std::count_if(edges.begin(), edges.begin() + mid, partially_local(seqs.front().id, seqs.back().id));

    return std::make_pair(true, "");
} // validate_edges

#endif // VALIDATE_EDGES_HPP
