/***
 *  $Id$
 **
 *  File: generate_edges.hpp
 *  Created: May 24, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2012-2015 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef GENERATE_EDGES_HPP
#define GENERATE_EDGES_HPP

#include <algorithm>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <vector>

#include <boost/container/vector.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <jaz/algorithm.hpp>
#include <mpix2/partition_balance.hpp>
#include <mpix2/sample_sort.hpp>
#include <mpix2/simple_partition.hpp>
#include <mpix2/write_cbuffer.hpp>

#include "../sequence_compact.hpp"

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"


#ifdef WITH_MPE
#include <mpix2/MPE_Log.hpp>
#endif // WITH_MPE


inline int partition_level(int size) {
    // ref is based on empirical data
    // depending on the number of processors 'size'
    // target blocks for balancing should be equal to
    // the largest partition divided by 'pos + 1'
    const unsigned int ref[] = { 6, 22, 84, 328, 1296, 5152, 20544, 82048, 327936, 1311232 };
    const int m = sizeof(ref);
    int pos = 0;
    for (; (pos < m) && (ref[pos] < size); ++pos);
    return (pos + 1);
} // partition_level


inline void extract_sketches(const AppConfig& opt, const SequenceList& SL, int i,
                             const std::vector<shingle_list_type>& shingles, std::vector<sketch_id>& sketch_list) {
    // get sketch-read pairs for given mod value
    unsigned int n = SL.seqs.size();

    for (unsigned int j = 0; j < n; ++j) {
        unsigned int id = SL.seqs[j].id;

        unsigned int l = shingles[j].size();
        unsigned int pos = sketch_list.size();

        for (unsigned int k = 0; k < l; ++k) {
            if (shingles[j][k] % opt.mod == i) {
                sketch_list.push_back(make_sketch_id(shingles[j][k], id, 0));
            }
        }

        l = sketch_list.size();
        for (unsigned int k = pos; k < l; ++k) sketch_list[k].size = (l - pos);
    } // for j
} // extract_sketches


inline unsigned int clean_sketches(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                                   std::vector<sketch_id>& sketch_list,
                                   std::vector<id_sketch>& rem_list) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // create rem_list and remove singletons
    if (rank == opt.dbg) report.stream << debug << "creating rem_list" << std::endl;

    std::sort(sketch_list.begin(), sketch_list.end());

    typedef std::vector<sketch_id>::iterator si_iterator;

    si_iterator iter = sketch_list.begin();
    si_iterator end = sketch_list.end();

    unsigned int max_part = 0;

    while (iter != end) {
        si_iterator temp = jaz::range(iter, end);
        if (temp - iter >= opt.cmax) {
            // add to the rem_list
            for (si_iterator i = iter; i != temp; ++i) {
                if (rem_list.empty() || (rem_list.back().id != i->id)) {
                    rem_list.push_back(make_id_sketch(i->id, i->sketch));
                }
            }

            // and remove from sketch_list
            std::copy(temp, end, iter);
            end -= (temp - iter);
        } else if (temp - iter == 1) {
            std::copy(temp, end, iter);
            end -= 1;
        } else {
            if (max_part < (temp - iter)) max_part = temp - iter;
            iter = temp;
        }
    } // while

    sketch_list.resize(end - sketch_list.begin());

    if (rank == opt.dbg) {
        report.stream << debug << "sketch_list: " << sketch_list.size()
                      << ", memory: " << to_size(sketch_list.size() * sizeof(sketch_id))
                      << "/" << to_size(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;

        report.stream << debug << "rem_list: " << rem_list.size()
                      << ", memory: " << to_size(rem_list.size() * sizeof(id_sketch))
                      << "/" << to_size(rem_list.capacity() * sizeof(id_sketch)) << std::endl;
    }

    return max_part;
} // clean_sketches


inline void balance_sketches(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                             unsigned int max_part, std::vector<sketch_id>& sketch_list) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // balance sketches
    if (rank == opt.dbg) report.stream << debug << "balancing sketches" << std::endl;

    MPI_Datatype MPI_SKETCH_ID;
    MPI_Type_contiguous(sizeof(sketch_id), MPI_BYTE, &MPI_SKETCH_ID);
    MPI_Type_commit(&MPI_SKETCH_ID);

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("decompose sketch_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    int part = 4 * opt.eps;

    if (opt.eps == 0) {
        // find optimal size of partitioning
        // proteins may require more complex strategy
        unsigned int gmax_part = 0;
        MPI_Allreduce(&max_part, &gmax_part, 1, MPI_UNSIGNED, MPI_MAX, comm);
        part = gmax_part / partition_level(size);
    }

    unsigned int part_id = rank * (std::numeric_limits<unsigned int>::max() / size) + 1;


    // std::ostringstream os;
    // print_sketch_tasks(os, sketch_list);
    // os << "\n";
    // mpix::write_cbuffer("sketch_list.balance0", os.str().c_str(), os.str().size(), comm);

    mpix::partition_balance(sketch_list, MPI_SKETCH_ID, comm);
    part_id = decompose_sketch_list(sketch_list, part, part_id);

    mpix::partition_balance(sketch_list, MPI_SKETCH_ID, comm);
    part_id = decompose_sketch_list(sketch_list, std::max(part / 4, 1024), part_id);

    // os.str("");
    // print_sketch_tasks(os, sketch_list);
    // os << "\n";
    // mpix::write_cbuffer("sketch_list.balance1", os.str().c_str(), os.str().size(), comm);

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

#ifdef WITH_MPE
    mpe_log.init("balance sketch_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    // here we actually do the balancing
    typedef std::vector<sketch_id>::iterator si_iterator;

    mpix::partition_balance(sketch_list, std::equal_to<sketch_id>(),
                            sketch_part_cost<si_iterator>, MPI_SKETCH_ID, 0, comm);

    // os.str("");
    // print_sketch_tasks(os, sketch_list);
    // os << "\n";
    // mpix::write_cbuffer("sketch_list.balance2", os.str().c_str(), os.str().size(), comm);

    if (rank == opt.dbg) {
        report.stream << debug << "sketch_list: " << sketch_list.size()
                      << ", memory: " << to_size(sketch_list.size() * sizeof(sketch_id))
                      << "/" << to_size(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
    }

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    MPI_Type_free(&MPI_SKETCH_ID);
} // balance_sketches


template <typename Sequence>
std::pair<bool, std::string> count_edges_support(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                                                 std::vector<sketch_id>& sketch_list,
                                                 Sequence& counts) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("get counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    // count required memory
    unsigned int count_alloc = 0;
    unsigned int bsz = 0;

    typedef std::vector<sketch_id>::iterator si_iterator;
    si_iterator iter = sketch_list.begin();

    // sketch list should be sorted/grouped by now
    while (iter != sketch_list.end()) {
        si_iterator temp = jaz::range(iter, sketch_list.end());
        int k = temp - iter;

        if (iter->sep == 0) bsz = nc2(k);
        else {
            if (iter->sep == k) bsz = k * k;
            else bsz = iter->sep * (k - iter->sep);
        }

        count_alloc += bsz;
        iter = temp;
    } // while


    // pre-allocate!
    if (rank == opt.dbg) {
        report.stream << debug << "allocating " << count_alloc << " counts" << std::endl;
    }

    try {
        // +1 just to be on the safe side :-)
        counts.reserve(count_alloc + 1);
    } catch (...) {
        return std::make_pair(false, "counts allocation failure,...");
    }

    if (rank == opt.dbg) {
        report.stream << debug << "counts, reserved: " << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }


    // generate counts
    if (rank == opt.dbg) {
        report.stream << debug << "generating counts" << std::endl;
    }

    iter = sketch_list.begin();

    while (iter != sketch_list.end()) {
        si_iterator temp = jaz::range(iter, sketch_list.end());

        if (iter->sep == 0) {
            // enumerate all pairs with the same sketch (triangle)
            for (si_iterator j = iter; j != temp - 1; ++j) {
                for (si_iterator k = j + 1; k != temp; ++k) {
                    if (j->id != k->id) counts.push_back(make_read_pair(*j, *k));
                } // for k
            } // for j
        } else {
            // enumerate all pairs with the same sketch (rectangle)
            si_iterator pos = iter + iter->sep;
            if (iter->sep == (temp - iter)) pos = iter;
            for (si_iterator j = iter; j != iter + iter->sep; ++j) {
                for (si_iterator k = pos; k != temp; ++k) {
                    if (j->id != k->id) counts.push_back(make_read_pair(*j, *k));
                } // for k
            } // for j
        } // if

        iter = temp;
    } // while

    return std::make_pair(true, "");
} // count_edges_support


template <typename Sequence>
std::pair<bool, std::string> compact_counts(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                                            const SequenceList& SL, const std::vector<read_pair>& edges,
                                            Sequence& counts) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("compact_counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    // read_xor almost uniformly divides edges between processors
    // matrix_block(SL.N, size) divides edges such that edges from
    // the same block in NxN matrix are assigned to the same processor
    std::transform(counts.begin(), counts.end(), counts.begin(), reads_xor);

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    if (rank == opt.dbg) report.stream << debug << "compact local" << std::endl;

    std::sort(counts.begin(), counts.end());
    counts.resize(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()) - counts.begin());

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    double t0 = 0;

    if (rank == opt.dbg) {
        report.stream << debug << "compact global" << std::endl;
        t0 = MPI_Wtime();
    }

    try {
        mpix::simple_partition(counts, hash_rp_score, MPI_READ_PAIR, comm);
        std::sort(counts.begin(), counts.end());

        counts.resize(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()) - counts.begin());
        if ((counts.size() > 0) && (counts.capacity() / counts.size() > 2)) counts.shrink_to_fit();
    } catch (...) {
        return std::make_pair(false, "counts compaction failure,...");
    }

    if (rank == opt.dbg) {
        double t1 = MPI_Wtime();
        report.stream << debug << "time: " << (t1 - t0) << "s" << std::endl;

        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    MPI_Type_free(&MPI_READ_PAIR);

    return std::make_pair(true, "");
} // compact_counts


template <typename Sequence>
std::pair<bool, std::string> compact_counts_lomem(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                                                  const SequenceList& SL, const std::vector<read_pair>& edges,
                                                  Sequence& counts) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("compact_counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    if (rank == opt.dbg) report.stream << debug << "compact local" << std::endl;

    std::sort(counts.begin(), counts.end());
    counts.resize(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()) - counts.begin());

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    double t0 = 0;

    if (rank == opt.dbg) {
        report.stream << debug << "compact global" << std::endl;
        t0 = MPI_Wtime();
    }

    try {
        boost::container::vector<read_pair> buf;
        boost::container::vector<read_pair> out;

        for (int i = 0; i < size; ++i) {
            if (rank == opt.dbg) report.stream << debug << "round: " << i << std::endl;

            tree_compact(counts, buf, hash_rp_xor, std::plus<read_pair>(), MPI_READ_PAIR, i, comm);
            if ((counts.size() > 0) && (counts.capacity() / counts.size() > 2)) counts.shrink_to_fit();

            if (i == rank) {
                out = buf;

                if (rank == opt.dbg) {
                    report.stream << debug << "out: " << out.size()
                                  << ", memory: " << to_size(out.size() * sizeof(read_pair))
                                  << "/" << to_size(out.capacity() * sizeof(read_pair)) << std::endl;
                }
            } // if i

            if (rank == opt.dbg) {
                report.stream << debug << "counts: " << counts.size()
                              << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                              << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
            }
        } // for i

        counts = boost::move(out);

    } catch (...) {
        return std::make_pair(false, "counts compaction failure,...");
    }

    if (rank == opt.dbg) {
        double t1 = MPI_Wtime();
        report.stream << debug << "time: " << (t1 - t0) << "s" << std::endl;

        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    MPI_Type_free(&MPI_READ_PAIR);

    return std::make_pair(true, "");
} // compact_counts_lomem


inline std::pair<bool, std::string> aggregate_aux_list(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                                                       std::vector<id_sketch>& rem_list) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("aggregate aux_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // sort rem_list for fast search
    MPI_Datatype MPI_ID_SKETCH;
    MPI_Type_contiguous(sizeof(id_sketch), MPI_BYTE, &MPI_ID_SKETCH);
    MPI_Type_commit(&MPI_ID_SKETCH);

    mpix::sample_sort(rem_list, MPI_ID_SKETCH, comm);

    // aggregate rem_list
    int rl_sz = rem_list.size();

    std::vector<int> rem_list_sz(size, 0);
    std::vector<int> rem_list_off(size, 0);

    MPI_Allgather(&rl_sz, 1, MPI_INT, &rem_list_sz[0], 1, MPI_INT, comm);

    unsigned int agg_rl_sz = std::accumulate(rem_list_sz.begin(), rem_list_sz.end(), 0);
    std::partial_sum(rem_list_sz.begin(), rem_list_sz.end() - 1, rem_list_off.begin() + 1);

    std::vector<id_sketch> agg_rem_list(agg_rl_sz);
    MPI_Allgatherv(&rem_list[0], rl_sz, MPI_ID_SKETCH, &agg_rem_list[0], &rem_list_sz[0], &rem_list_off[0], MPI_ID_SKETCH, comm);

    if (rank == opt.dbg) {
        report.stream << debug << "agg_rem_list: " << agg_rl_sz
                      << ", memory: " << to_size(agg_rl_sz * sizeof(id_sketch))
                      << "/" << to_size(agg_rem_list.capacity() * sizeof(id_sketch)) << std::endl;
    }

    MPI_Type_free(&MPI_ID_SKETCH);

    rem_list = agg_rem_list;

    return std::make_pair(true, "");
}; // aggregate_aux_list


template <typename Sequence>
void update_counts(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                   Sequence& counts, const std::vector<id_sketch>& rem_list) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("update counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == opt.dbg) report.stream << debug << "updating counts" << std::endl;

    typename Sequence::iterator first = counts.begin();
    typename Sequence::iterator last = counts.end();

    typedef std::vector<id_sketch>::const_iterator vp_iterator;

    std::pair<vp_iterator, vp_iterator> r0;
    std::pair<vp_iterator, vp_iterator> r1;

    read_pair rp;
    double jmin = opt.jmin;

    for (; first != last; ++first) {
        rp = kmer_fraction(*first);
        if (rp.count > jmin) continue;

        r0 = std::equal_range(rem_list.begin(), rem_list.end(), make_id_sketch(first->id0, 0), id_compare2);
        r1 = std::equal_range(r0.second, rem_list.end(), make_id_sketch(first->id1, 0), id_compare2);

        first->count += jaz::intersection_size(r0.first, r0.second, r1.first, r1.second, sketch_compare2);
    }
} // update_counts


template <typename Sequence>
void filter_counts(const AppConfig& opt, Reporter& report, MPI_Comm comm, Sequence& counts) {
#ifdef WITH_MPE
    mpe_log.init("filter counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // approximate kmer fraction and filter edges
    // count is now approximated kmer fraction
    if (rank == opt.dbg) report.stream << debug << "filtering counts" << std::endl;

    std::transform(counts.begin(), counts.end(), counts.begin(), kmer_fraction);
    counts.resize(std::remove_if(counts.begin(), counts.end(), not_similar(opt.jmin)) - counts.begin());

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << to_size(counts.size() * sizeof(read_pair))
                      << "/" << to_size(counts.capacity() * sizeof(read_pair)) << std::endl;
    }
} // filter_counts


template <typename Sequence>
std::pair<bool, std::string> merge_edges(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                                         Sequence& counts, std::vector<read_pair>& edges) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("merge edges", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == opt.dbg) report.stream << debug << "merging edges" << std::endl;

    unsigned int edges_sz = edges.size();
    unsigned int counts_sz = counts.size();

    std::sort(counts.begin(), counts.end());

    try {
        if (edges.empty()) {
            edges.resize(counts.size());
            std::copy(counts.begin(), counts.end(), edges.begin());
        } else {
            int i = 0;
            int j = 0;

            while ((i < counts_sz) && (j != edges_sz)) {
                if (counts[i] < edges[j]) edges.push_back(counts[i++]);
                else if (edges[j] < counts[i]) j++;
                else { i++; j++; }
            }

            std::copy(counts.begin() + i, counts.end(), std::back_inserter(edges));
            std::inplace_merge(edges.begin(), edges.begin() + edges_sz, edges.end());
        }
    } catch (...) {
        return std::make_pair(false, "edge merging failure,...");
    }

    if (rank == opt.dbg) {
        report.stream << debug << "edges: " << edges.size()
                      << ", +" << (edges.size() - edges_sz)
                      << " memory: " << to_size(edges.size() * sizeof(read_pair))
                      << "/" << to_size(edges.capacity() * sizeof(read_pair)) << std::endl;
    }

    return std::make_pair(true, "");
} // merge_edges


inline std::pair<bool, std::string> extract_edges(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
                                                  const SequenceList& SL, std::vector<sketch_id>& sketch_list,
                                                  std::vector<read_pair>& edges) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    bool res = false;
    std::string err = "";

    boost::container::vector<read_pair> counts;
    std::vector<id_sketch> rem_list;

    // remove singletons and frequent sketches
    unsigned int max_part = clean_sketches(opt, report, comm, sketch_list, rem_list);

    // balance sketches such that the total number of contributing pairs
    // generated per processor is more or less the same on every rank
    balance_sketches(opt, report, comm, max_part, sketch_list);

    // sketch_list is now sorted!!!

    // enumerate all pairs that will add to score for edge
    boost::tie(res, err) = count_edges_support(opt, report, comm, sketch_list, counts);
    if (res == false) return std::make_pair(false, err);

    sketch_list.clear();

    // compact counts per edge
    boost::tie(res, err) = compact_counts(opt, report, comm, SL, edges, counts);
    if (res == false) return std::make_pair(false, err);

    // prepare auxiliary list
    boost::tie(res, err) = aggregate_aux_list(opt, report, comm, rem_list);
    if (res == false) return std::make_pair(false, err);

    // update the final count per edge
    update_counts(opt, report, comm, counts, rem_list);

    // removed pairs without sufficient support
    filter_counts(opt, report, comm, counts);
    counts.shrink_to_fit();

    // add to candidate edges
    boost::tie(res, err) = merge_edges(opt, report, comm, counts, edges);
    if (res == false) return std::make_pair(false, err);

    return std::make_pair(true, "");
} // extract_edges


inline std::pair<bool, std::string> generate_edges(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
                                                   const SequenceList& SL, const std::vector<shingle_list_type>& shingles,
                                                   std::vector<read_pair>& edges) {
    report << step << "extracting candidate edges..." << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    std::vector<sketch_id> sketch_list;

    unsigned int n = SL.seqs.size();
    unsigned int end = std::min(opt.iter, opt.mod);

    // register MPI type for sketch_id pair
    MPI_Datatype MPI_SKETCH_ID;
    MPI_Type_contiguous(sizeof(sketch_id), MPI_BYTE, &MPI_SKETCH_ID);
    MPI_Type_commit(&MPI_SKETCH_ID);

    report << info << "running " << end << " iteration(s)..." << std::endl;

    bool res = false;
    std::string err = "";

    // main loop (we look for several evidences to get edge)
    for (unsigned int i = 0; i < end; ++i) {
        if (opt.dbg < 0) report << "." << std::flush;
        if (rank == opt.dbg) report.stream << debug << "iteration: " << i << std::endl;

#ifdef WITH_MPE
        mpix::MPE_Log mpe_log("generate sketch_list", "red");
        mpe_log.start();
#endif // WITH_MPE

        extract_sketches(opt, SL, i, shingles, sketch_list);

        if (rank == opt.dbg) {
            report.stream << debug << "sketch_list: " << sketch_list.size()
                          << ", memory: " << to_size(sketch_list.size() * sizeof(sketch_id))
                          << "/" << to_size(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
        }

#ifdef WITH_MPE
        mpe_log.stop();
#endif // WITH_MPE


#ifdef WITH_MPE
        mpe_log.init("group sketch_list", "red");
        mpe_log.start();
#endif // WITH_MPE

        // group globally sketch list
        if (rank == opt.dbg) report.stream << debug << "grouping sketches" << std::endl;

        mpix::simple_partition(sketch_list, hash_sketch_id2, MPI_SKETCH_ID, comm);

        if (rank == opt.dbg) {
            report.stream << debug << "sketch_list: " << sketch_list.size()
                          << ", memory: " << to_size(sketch_list.size() * sizeof(sketch_id))
                          << "/" << to_size(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
        }

#ifdef WITH_MPE
        mpe_log.stop();
#endif // WITH_MPE


        // add to edges read pairs with common sketches
        // this is the key ingredient
        unsigned int edges_sz = edges.size();

        boost::tie(res, err) = extract_edges(opt, log, report, comm, SL, sketch_list, edges);
        if (res == false) return std::make_pair(false, err);

        if (opt.dbg < 0) report << ":" << std::flush;

        // check if we are getting new edges
        unsigned int ed = edges.size() - edges_sz;
        MPI_Allreduce(&ed, &edges_sz, 1, MPI_UNSIGNED, MPI_SUM, comm);

        if (rank == opt.dbg) {
            report.stream << debug << "total edges added: " << edges_sz
                          << " memory: " << to_size(edges_sz * sizeof(read_pair))
                          << std::endl;
            report.stream << debug << std::endl;
        }
    } // for i

    MPI_Type_free(&MPI_SKETCH_ID);

    if (opt.dbg < 0) report << std::endl;

    return std::make_pair(true, "");
} // generate_edges

#endif // GENERATE_EDGES_HPP
