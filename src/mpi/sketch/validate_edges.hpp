/***
 *  $Id$
 **
 *  File: validate_edges.hpp
 *  Created: May 29, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2012-2014 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef VALIDATE_EDGES_HPP
#define VALIDATE_EDGES_HPP

#include <string>
#include <vector>

#include <mpix2/simple_partition.hpp>

#include "../StaticWSQueue.hpp"

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "create_smatrix.hpp"

#include "iomanip.hpp"


#ifdef WITH_MPE
#include <mpix2/MPE_Log.hpp>
#endif // WITH_MPE


class compare_method {
public:
    explicit compare_method(const AppConfig& opt) : method_(opt.method) {
        if (opt.method == 0) {
            kf_ = bio::kmer_fraction(opt.kmer, opt.is_dna);
        } else {
            bio::scoring_matrix sm;
            int g, h;

            create_smatrix(opt.gaps, opt.is_dna, sm, g, h);

            if ((opt.method == 1) || (opt.method == 3)) {
                cfe_align_ = bio::free_global_alignment(sm, g, h);
            } else if ((opt.method == 2) || (opt.method == 4)) {
                align_ = bio::global_alignment(sm, g, h);
            } else if (opt.method == 5) loc_align_ = bio::local_alignment(sm, g, h);
            else balign_ = bio::banded_global_alignment(sm, g, h, 3 * opt.kmer);
        }
    } // compare_method

    boost::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
        const std::string* sa = &s0;
        const std::string* sb = &s1;

        if (s1.size() < s0.size()) std::swap(sa, sb);

        boost::tuple<int, int, int> res = boost::make_tuple(-1, -1, -1);

        if (method_ == 0) res = kf_(*sa, *sb);
        else if ((method_ == 1) || (method_ == 3)) res = cfe_align_(*sa, *sb);
        else if ((method_ == 2) || (method_ == 4)) res = res = align_(*sa, *sb);
        else if (method_ == 5) res = loc_align_(*sa, *sb);
        else res = balign_(*sa, *sb);

        // correction to get score for CD-HIT identity score
        if ((method_ == 1) || (method_ == 2) || (method_ == 6)) boost::get<1>(res) = std::min(s0.size(), s1.size());

        return res;
    } // operator()

private:
    bio::kmer_fraction kf_;
    bio::global_alignment align_;
    bio::free_global_alignment cfe_align_;
    bio::local_alignment loc_align_;
    bio::banded_global_alignment balign_;

    unsigned int method_;

}; // struct compare_method


inline bool block_compare(const std::pair<unsigned int, unsigned int>& p1, const std::pair<unsigned int, unsigned int>& p2) {
    return (p2.second - p2.first) < (p1.second - p1.first);
} // block_compare


inline std::pair<bool, std::string> validate_edges(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
                                                   const SequenceList& SL, SequenceRMA& rma_seq,
                                                   std::vector<read_pair>& edges) {
    report << step << "validating edges:" << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // rebalance graph
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    mpix::simple_partition(edges, hash_rp_id0or1(SL.N, size), MPI_READ_PAIR, comm);

    MPI_Type_free(&MPI_READ_PAIR);

    // divide into local and remote
    unsigned int mid =
        std::partition(edges.begin(), edges.end(), local(SL.seqs.front().id, SL.seqs.back().id)) - edges.begin();

    // local edges are easy :-)
    report << info << "processing local edges..." << std::endl;

    sequence_compare<compare_method> ident(SL, SL.seqs.front().id, compare_method(opt));
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

    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar2(opt.method, opt.level)), edges.end());
    if (opt.factor == false) std::transform(edges.begin(), edges.end(), edges.begin(), read_pair_count(opt.method));

    return std::make_pair(true, "");
} // validate_edges


inline std::pair<bool, std::string> validate_edges_ws(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
                                                      const SequenceList& SL, SequenceRMA& rma_seq,
                                                      std::vector<read_pair>& edges) {
    report << step << "validating edges:" << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // rebalance graph
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    mpix::simple_partition(edges, hash_rp_id0or1(SL.N, size), MPI_READ_PAIR, comm);

    MPI_Type_free(&MPI_READ_PAIR);

    report << info << "creating tasks..." << std::endl;

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("transform edges", "green");
    mpe_log.start();
#endif // WITH_MPE

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

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

#ifdef WITH_MPE
    mpe_log.init("prepare task queue", "green");
    mpe_log.start();
#endif // WITH_MPE

    // get tasks list
    typedef StaticWSQueue<read_pair> ws_queue_type;
    const unsigned int SBLOCK = 1024;

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

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    // process local edges
    report << info << "processing local edges..." << std::endl;

#ifdef WITH_MPE
    mpe_log.init("process local edges", "green");
    mpe_log.start();
#endif // WITH_MPE

    sequence_compare<compare_method> ident(SL, SL.seqs.front().id, compare_method(opt));
    edges.resize(mid);

    std::transform(edges.begin(), edges.end(), edges.begin(), ident);
    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar2(opt.method, opt.level)), edges.end());

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    // process tasks
    report << info << "processing remaining edges, be patient..." << std::endl;

#ifdef WITH_MPE
    mpe_log.init("process semi-local edges", "green");
    mpe_log.start();
#endif // WITH_MPE

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
            wsq.progress();

            const std::string& s = sv[rank_id[j].second - start_id];
            const read_pair& e = first[j];

            if (i2r(e.id0) == rank) edges.push_back(ident(e, s));
            else edges.push_back(ident(s, e));
        } // for j

        wsq.progress();
    } // while wsq.get

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    // we start stealing process
    report << info << "stealing work..." << std::endl;

#ifdef WITH_MPE
    mpe_log.init("process remote edges", "green");
    mpe_log.start();
#endif // WITH_MPE

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

            wsq.progress();
        } // for

        delete[] sfirst;
    } // while wsq.steal

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

#ifdef WITH_MPE
    mpe_log.init("final cleaning", "green");
    mpe_log.start();
#endif // WITH_MPE

    edges.erase(std::remove_if(edges.begin(), edges.end(), not_similar2(opt.method, opt.level)), edges.end());
    if (opt.factor == false) std::transform(edges.begin(), edges.end(), edges.begin(), read_pair_count(opt.method));

    wsq.finalize();

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    return std::make_pair(true, "");
} // validate_edges_ws

#endif // VALIDATE_EDGES_HPP
