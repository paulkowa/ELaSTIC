/***
 *  $Id$
 **
 *  File: generate_edges.hpp
 *  Created: May 24, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2012-2014 Jaroslaw Zola
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
#include <sstream>
#include <string>
#include <vector>

#include <boost/container/vector.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <jaz/algorithm.hpp>
#include <mpix2/partition_balance.hpp>
#include <mpix2/sample_sort.hpp>
#include <mpix2/simple_partition.hpp>

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"


#ifdef WITH_MPE
#include <mpix2/MPE_Log.hpp>
#endif // WITH_MPE


template <typename Int> inline Int nc2(Int n) { return (n * (n - 1)) >> 1; }


inline std::string size_str(unsigned int sz) {
    std::ostringstream ss;
    if (sz < 1024) ss << sz << "B";
    else if (sz < (1024 * 1024)) {
        ss << (static_cast<double>(sz) / 1024) << "KB";
    } else {
        ss << (static_cast<double>(sz) / (1024 * 1024)) << "MB";
    }
    return ss.str();
} // size_str


inline int partition_level(int size) {
    const unsigned int ref[] = { 6, 22, 84, 328, 1296, 5152, 20544, 82048, 327936, 1311232 };
    int pos = 0;
    for (; (pos < sizeof(ref)) && (ref[pos] < size); ++pos);
    return (pos + 1);
} // partition_level


template <typename Iter>
void update_counts(Iter first, Iter last, const std::vector<id_sketch>& rem_list, double jmin) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("update_counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    typedef std::vector<id_sketch>::const_iterator vp_iterator;

    std::pair<vp_iterator, vp_iterator> r0;
    std::pair<vp_iterator, vp_iterator> r1;

    read_pair rp;

    for (; first != last; ++first) {
        rp = kmer_fraction(*first);
        if (rp.count > jmin) continue;

        r0 = std::equal_range(rem_list.begin(), rem_list.end(), make_id_sketch(first->id0, 0), id_compare2);
        r1 = std::equal_range(r0.second, rem_list.end(), make_id_sketch(first->id1, 0), id_compare2);

        first->count += jaz::intersection_size(r0.first, r0.second, r1.first, r1.second, sketch_compare2);
    }
} // update_counts


class in_edges {
public:
    in_edges(const std::vector<read_pair>& edges) : edges_(edges) { }

    bool operator()(const read_pair& rp) const {
        return std::binary_search(edges_.begin(), edges_.end(), rp);
    } // operator()

private:
    const std::vector<read_pair>& edges_;
}; // class in_edges


// if compact_counts changes the edge merging step must be updated!!!
template <typename Sequence>
void compact_counts(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                    unsigned int N, const std::vector<read_pair>& edges, Sequence& counts) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("compact_counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    std::transform(counts.begin(), counts.end(), counts.begin(), matrix_block(N, size));

    Sequence data;
    data.reserve(counts.size() + (counts.size() >> 2));

    const int DS = 10;
    unsigned int bsz = std::ceil(1.0 * counts.size() / DS);

    std::vector<read_pair> buf;
    typename Sequence::iterator iter;

    if (rank == opt.dbg) {
        report.stream << debug << "DS: " << DS << " bsz: " << bsz << std::endl;
    }

    // this will work as long as counts are fairly randomly distributed
    for (int i = 0; i < DS; ++i) {
        iter = std::max(counts.begin(), counts.end() - bsz);
        buf.resize(counts.end() - iter);

        std::copy(iter, counts.end(), buf.begin());
        counts.resize(iter - counts.begin());

        mpix::simple_partition(buf, hash_read_pair3, MPI_READ_PAIR, comm);
        buf.resize(std::remove_if(buf.begin(), buf.end(), in_edges(edges)) - buf.begin());
        std::copy(buf.begin(), buf.end(), std::back_inserter(data));

        if (rank == opt.dbg) {
            unsigned int mem = (buf.size() + counts.size() + data.size()) * sizeof(read_pair);
            unsigned int mem_res = (buf.capacity() + counts.capacity() + data.capacity()) * sizeof(read_pair);

            report.stream << debug << i
                          << " counts: " << counts.size() << "/" << counts.capacity()
                          << "\tbuf: " << buf.size() << "/" << buf.capacity()
                          << "\tdata: " << data.size() << "/" << data.capacity()
                          << "\tmemory: " << size_str(mem) << "/" << size_str(mem_res) << std::endl;
        }
    } // for i

    counts = boost::move(data);

    // run actual compaction
    if (rank == opt.dbg) report.stream << debug << "compact" << std::endl;

    std::sort(counts.begin(), counts.end());
    counts.resize(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()) - counts.begin());

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << size_str(counts.size() * sizeof(read_pair))
                      << "/" << size_str(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    MPI_Type_free(&MPI_READ_PAIR);

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE
} // compact_counts


template <typename Sequence>
void aggregate_list(const AppConfig& opt, Reporter& report, MPI_Comm comm,
                    const Sequence& counts, std::vector<id_sketch>& rem_list) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("aggregate_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    MPI_Datatype MPI_ID_SKETCH;
    MPI_Type_contiguous(sizeof(id_sketch), MPI_BYTE, &MPI_ID_SKETCH);
    MPI_Type_commit(&MPI_ID_SKETCH);

    // sort rem_list for fast search
    if (rank == opt.dbg) report.stream << debug << "boundaries" << std::endl;

    mpix::sample_sort(rem_list, MPI_ID_SKETCH, comm);

    // distribute current bin boundaries
    std::vector<unsigned int> ranges(size, 0);

    int lid = rem_list.empty() ? 0 : rem_list.back().id;
    MPI_Allgather(&lid, 1, MPI_UNSIGNED, &ranges[0], 1, MPI_UNSIGNED, comm);

    // check if global rem_list is not empty
    if (std::accumulate(ranges.begin(), ranges.end(), 0) == 0) return;

    // rebalance by id overlapping lists
    unsigned int sexchange = 0;
    unsigned int rexchange = 0;

    if ((rank > 0) && (!rem_list.empty())) {
        id_sketch is = make_id_sketch(rem_list.front().id, 0);

        if (ranges[rank - 1] == rem_list.front().id) {
            sexchange = std::upper_bound(rem_list.begin(), rem_list.end(), is, id_compare2) - rem_list.begin();
        }
    }

    const int EXCHANGE_TAG = 111;
    MPI_Status stat;

    unsigned int rl_sz = rem_list.size();

    if (rank == 0) {
        MPI_Recv(&rexchange, 1, MPI_UNSIGNED, rank + 1, EXCHANGE_TAG, comm, &stat);
        rem_list.resize(rl_sz + rexchange);
        MPI_Recv(&rem_list[rl_sz], rexchange, MPI_ID_SKETCH, rank + 1, EXCHANGE_TAG, comm, &stat);
    } else if (rank == size - 1) {
        MPI_Send(&sexchange, 1, MPI_UNSIGNED, rank - 1, EXCHANGE_TAG, comm);
        MPI_Send(&rem_list[0], sexchange, MPI_ID_SKETCH, rank - 1, EXCHANGE_TAG, comm);
        rem_list.erase(rem_list.begin(), rem_list.begin() + sexchange);
    } else {
        MPI_Sendrecv(&sexchange, 1, MPI_UNSIGNED, rank - 1, EXCHANGE_TAG,
                     &rexchange, 1, MPI_UNSIGNED, rank + 1, EXCHANGE_TAG, comm, &stat);
        rem_list.resize(rl_sz + rexchange);
        MPI_Sendrecv(&rem_list[0], sexchange, MPI_ID_SKETCH, rank - 1, EXCHANGE_TAG,
                     &rem_list[rl_sz], rexchange, MPI_ID_SKETCH, rank + 1, EXCHANGE_TAG, comm, &stat);
        rem_list.erase(rem_list.begin(), rem_list.begin() + sexchange);
    }

    // now we can send info about new limits and sizes
    // we add 1 to get one-past-end element
    std::vector<unsigned int> rsizes(size, 0);

    lid = rem_list.empty() ? 0 : rem_list.back().id + 1;
    MPI_Allgather(&lid, 1, MPI_UNSIGNED, &ranges[0], 1, MPI_UNSIGNED, comm);

    lid = rem_list.size();
    MPI_Allgather(&lid, 1, MPI_UNSIGNED, &rsizes[0], 1, MPI_UNSIGNED, comm);

    // correct for empty range
    for (int i = 1; i < size; ++i) if (ranges[i] == 0) ranges[i] = ranges[i - 1];

    // find which bins we need and their size
    std::vector<unsigned int>::iterator iter;

    std::vector<int> rbins(size, -1);
    std::vector<int> dbins(size, -1);

    for (int i = 0; i < counts.size(); ++i) {
        iter = std::lower_bound(ranges.begin(), ranges.end(), counts[i].id0);
        rbins[iter - ranges.begin()] = rank;
        iter = std::lower_bound(iter, ranges.end(), counts[i].id1);
        rbins[iter - ranges.begin()] = rank;
    }

    if (rank == opt.dbg) report.stream << debug << "aggregation" << std::endl;

    unsigned int rb = 0;
    unsigned int S = 0;

    for (int i = 0; i < size; ++i) if (rbins[i] == rank) { S += rsizes[i]; rb++; }

    // tell others
    MPI_Alltoall(&rbins[0], 1, MPI_INT, &dbins[0], 1, MPI_INT, comm);

    // we have to send our chunk to dbins[] != -1
    std::vector<id_sketch> agg_rem_list(S);

    if (rank == opt.dbg) {
        report.stream << debug << "agg_rem_list: " << S
                      << ", memory: " << size_str(S * sizeof(id_sketch))
                      << "/" << size_str(agg_rem_list.capacity() * sizeof(id_sketch))
                      << " parts: " << rb << std::endl;
    }

    std::vector<int> scounts(size, 0);
    std::vector<int> sdispl(size, 0);

    for (int i = 0; i < size; ++i) if (dbins[i] != -1) scounts[i] = rem_list.size();

    std::vector<int> rcounts(size, 0);
    std::vector<int> rdispl(size, 0);

    for (int i = 0; i < size; ++i) if (rbins[i] == rank) rcounts[i] = rsizes[i];
    std::partial_sum(rcounts.begin(), rcounts.end() - 1, rdispl.begin() + 1);

    MPI_Alltoallv(&rem_list[0], &scounts[0], &sdispl[0], MPI_ID_SKETCH,
                  &agg_rem_list[0], &rcounts[0], &rdispl[0], MPI_ID_SKETCH, comm);

    MPI_Type_free(&MPI_ID_SKETCH);

    rem_list = agg_rem_list;
}; // aggregate_list


inline std::pair<bool, std::string> extract_pairs(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
                                                  std::vector<sketch_id>& sketch_list, unsigned int N,
                                                  std::vector<read_pair>& edges) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // create rem_list and remove singletons
    if (rank == opt.dbg) report.stream << debug << "creating rem_list" << std::endl;

    std::sort(sketch_list.begin(), sketch_list.end());
    std::vector<id_sketch> rem_list;

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
        report.stream << debug << "rem_list: " << rem_list.size()
                      << ", memory: " << size_str(rem_list.size() * sizeof(id_sketch))
                      << "/" << size_str(rem_list.capacity() * sizeof(id_sketch)) << std::endl;

        report.stream << debug << "sketch_list: " << sketch_list.size()
                      << ", memory: " << size_str(sketch_list.size() * sizeof(sketch_id))
                      << "/" << size_str(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
    }

    // balance sketches
    if (rank == opt.dbg) report.stream << debug << "balancing sketches" << std::endl;

    MPI_Datatype MPI_SKETCH_ID;
    MPI_Type_contiguous(sizeof(sketch_id), MPI_BYTE, &MPI_SKETCH_ID);
    MPI_Type_commit(&MPI_SKETCH_ID);

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("decompose sketch_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    int part = 2 * opt.eps;

    if (opt.eps == 0) {
        // find optimal size of partitioning
        // sizes based on empirical data
        unsigned int gmax_part = 0;
        MPI_Allreduce(&max_part, &gmax_part, 1, MPI_UNSIGNED, MPI_MAX, comm);
        part = gmax_part / partition_level(size);
    }

    unsigned int part_id = rank * (std::numeric_limits<unsigned int>::max() / size) + 1;

    mpix::partition_balance(sketch_list, MPI_SKETCH_ID, comm);
    part_id = decompose_sketch_list(sketch_list, part, part_id);

    mpix::partition_balance(sketch_list, MPI_SKETCH_ID, comm);
    part_id = decompose_sketch_list(sketch_list, std::max(part / 4, 500), part_id);

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

#ifdef WITH_MPE
    mpe_log.init("balance sketch_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    mpix::partition_balance(sketch_list, std::equal_to<sketch_id>(),
                            sketch_part_cost<si_iterator>, MPI_SKETCH_ID, 0, comm);

    if (rank == opt.dbg) {
        report.stream << debug << "sketch_list: " << sketch_list.size()
                      << ", memory: " << size_str(sketch_list.size() * sizeof(sketch_id))
                      << "/" << size_str(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
    }

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    MPI_Type_free(&MPI_SKETCH_ID);

#ifdef WITH_MPE
    mpe_log.init("get counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    // get local counts
    boost::container::vector<read_pair> counts;

    // count required memory
    unsigned int count_alloc = 0;
    unsigned int bsz = 0;

    iter = sketch_list.begin();

    while (iter != sketch_list.end()) {
        si_iterator temp = jaz::range(iter, sketch_list.end());

        if (iter->sep == 0) bsz = nc2(temp - iter);
        else {
            if (iter->sep == (temp - iter)) bsz = (temp - iter) * (temp - iter);
            else bsz = iter->sep * ((temp - iter) - iter->sep);
        }

        count_alloc += bsz;
        iter = temp;
    } // while

    if (rank == opt.dbg) {
        report.stream << debug << "allocating " << count_alloc << " counts" << std::endl;
    }

    try {
        counts.reserve(count_alloc);
    } catch (...) {
        return std::make_pair(false, "allocation failure,...");
    }

    if (rank == opt.dbg) {
        report.stream << debug << "counts, reserved: " << counts.capacity() * sizeof(read_pair) << std::endl;
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

    sketch_list.clear();

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    // perform counts compaction
    if (rank == opt.dbg) report.stream << debug << "compacting counts" << std::endl;

    try {
        compact_counts(opt, report, comm, N, edges, counts);
    } catch (...) {
        return std::make_pair(false, "compaction failure,...");
    }

    // aggregate rem_list
    if (rank == opt.dbg) report.stream << debug << "aggregating rem_list" << std::endl;

    try {
        aggregate_list(opt, report, comm, counts, rem_list);
    } catch (...) {
        return std::make_pair(false, "aggregation failure,...");
    }

    // update counts
    if (rank == opt.dbg) report.stream << debug << "updating counts" << std::endl;

    update_counts(counts.begin(), counts.end(), rem_list, opt.jmin);

#ifdef WITH_MPE
    mpe_log.init("filter candidate edges", "red");
    mpe_log.start();
#endif // WITH_MPE

    // approximate kmer fraction and filter edges
    // count is now approximated kmer fraction
    if (rank == opt.dbg) report.stream << debug << "filtering edges" << std::endl;

    std::transform(counts.begin(), counts.end(), counts.begin(), kmer_fraction);
    counts.resize(std::remove_if(counts.begin(), counts.end(), not_similar(opt.jmin)) - counts.begin());

    if (rank == opt.dbg) {
        report.stream << debug << "counts: " << counts.size()
                      << ", memory: " << size_str(counts.size() * sizeof(read_pair))
                      << "/" << size_str(counts.capacity() * sizeof(read_pair)) << std::endl;
    }

    // add new edges
    if (rank == opt.dbg) report.stream << debug << "merging edges" << std::endl;

    // if compact_counts changes this has to be updated!!!
    std::copy(counts.begin(), counts.end(), std::back_inserter(edges));
    std::sort(edges.begin(), edges.end());

    if (rank == opt.dbg) {
        report.stream << debug << "edges: " << edges.size()
                      << ", +" << counts.size()
                      << " memory: " << size_str(edges.size() * sizeof(read_pair))
                      << "/" << size_str(edges.capacity() * sizeof(read_pair)) << std::endl;
    }

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    return std::make_pair(true, "");
} // extract_pairs


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

    // number of edge improvements
    unsigned int tries = 0;

    // main loop (we look for several evidences to get edge)
    for (unsigned int i = 0; i < end; ++i) {
        if (opt.dbg < 0) report << "." << std::flush;
        if (rank == opt.dbg) report.stream << debug << "iteration: " << i << std::endl;

#ifdef WITH_MPE
        mpix::MPE_Log mpe_log("generate sketch_list", "red");
        mpe_log.start();
#endif // WITH_MPE

        // get sketch-read pairs for given mod value
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

        if (rank == opt.dbg) {
            report.stream << debug << "sketch_list: " << sketch_list.size()
                          << ", memory: " << size_str(sketch_list.size() * sizeof(sketch_id))
                          << "/" << size_str(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
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
                          << ", memory: " << size_str(sketch_list.size() * sizeof(sketch_id))
                          << "/" << size_str(sketch_list.capacity() * sizeof(sketch_id)) << std::endl;
        }

#ifdef WITH_MPE
        mpe_log.stop();
#endif // WITH_MPE

        /* we do re-balancing anyway
        if (sketch_list.empty()) {
            report.critical << warning << "{" << rank << "}" << " empty list of sketches!" << std::endl;
        }
        */

        // add to edges read pairs with common sketches
        unsigned int edges_sz = edges.size();

        boost::tie(res, err) = extract_pairs(opt, log, report, comm, sketch_list, SL.N, edges);
        if (res == false) return std::make_pair(false, err);

        if (opt.dbg < 0) report << ":" << std::flush;

        // check if we are getting new edges
        unsigned int ed = edges.size() - edges_sz;
        MPI_Allreduce(&ed, &edges_sz, 1, MPI_UNSIGNED, MPI_SUM, comm);

        if (edges_sz == 0) tries++; else tries = 0;
        if (tries == 2) break;
    } // for i

    MPI_Type_free(&MPI_SKETCH_ID);

    if (opt.dbg < 0) report << std::endl;

    return std::make_pair(true, "");
} // generate_edges

#endif // GENERATE_EDGES_HPP
