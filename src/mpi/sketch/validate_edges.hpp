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

    std::vector<read_pair>::iterator mid, last;

    last = std::partition(edges.begin(), edges.end(), not_local(seqs.front().id, seqs.back().id));
    mid = std::partition(edges.begin(), last, not_collocated(log.input, size));

    unsigned int loc[2];
    loc[0] = edges.end() - last;
    loc[1] = last - mid;

    unsigned int min[2] = { 0, 0 };
    unsigned int max[2] = { 0, 0 };
    unsigned int tot[2] = { 0, 0 };

    MPI_Reduce(loc, min, 2, MPI_UNSIGNED, MPI_MIN, 0, comm);
    MPI_Reduce(loc, max, 2, MPI_UNSIGNED, MPI_MAX, 0, comm);

    report << info << "local edges distribution [" << min[0] << "," << max[0] << "]" << std::endl;
    report << info << "remote edges distribution [" << min[1] << "," << max[1] << "]" << std::endl;

    // first local edges (i.e. reads and edge are available)
    std::transform(last, edges.end(), last, identity(seqs, opt.method, opt.kmer));
    edges.erase(std::remove_if(last, edges.end(), not_similar(opt.level)), edges.end());

    report << info << "local edges processed" << std::endl;

    // now we group edges on the processors that can process them
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    read_pair* efirst = 0;
    read_pair* elast = 0;

    boost::tie(efirst, elast) = mpix::data_bucketing(mid, last, read_pair_local_hash(log.input, size),
						     MPI_READ_PAIR, 0, comm);
    edges.erase(mid, last);

    // we process received batch
    std::transform(efirst, elast, efirst, identity(seqs, opt.method, opt.kmer));
    elast = std::remove_if(efirst, elast, not_similar(opt.level));

    std::copy(efirst, elast, std::back_inserter(edges));

    delete[] efirst;

    report << info << "remote edges processed" << std::endl;

    MPI_Type_free(&MPI_READ_PAIR);

    return std::make_pair(true, "");
} // validate_edges

#endif // VALIDATE_EDGES_HPP
