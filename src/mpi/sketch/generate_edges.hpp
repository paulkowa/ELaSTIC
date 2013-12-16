/***
 *  $Id$
 **
 *  File: generate_edges.hpp
 *  Created: May 24, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the License.
 *  See accompanying LICENSE.
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

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <jaz/algorithm.hpp>
#include <mpix2/partition_balance.hpp>
#include <mpix2/simple_partition.hpp>

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"

#ifdef WITH_MPE
#include <mpix2/MPE_Log.hpp>
#endif // WITH_MPE


inline int partition_level(int size) {
    const unsigned int ref[] = { 6, 22, 84, 328, 1296, 5152, 20544, 82048, 327936, 1311232 };
    int pos = 0;
    for (; (pos < sizeof(ref)) && (ref[pos] < size); ++pos);
    return (pos + 1);
} // partition_level


template <typename T, typename U>
inline bool less1st(const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
    return (p1.first < p2.first);
} // less1st

template <typename T, typename U>
inline bool less2nd(const std::pair<T, U>& p1, const std::pair<T, U>& p2) {
    return (p1.second < p2.second);
} // less2nd


template <typename Iter>
inline void update_counts2(Iter first, Iter last, std::vector<int>& rem_list, double jmin) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("update_counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    std::vector<int> rem_index(1, 0);
    for (int i = 0; i != rem_list.size(); ++i) if (rem_list[i] == -1) rem_index.push_back(i);

    std::vector<std::pair<int, int> > lst;

    int l = rem_index.size() - 1;
    read_pair rp;

    for (int i = l; i != 0; --i) {
	int pos = rem_index[i - 1];
	int end = rem_index[i];
	for (int j = pos; j != end; ++j) if (rem_list[j] != -1) lst.push_back(std::make_pair(rem_list[j], i));
	rem_list.resize(pos);
    } // for i

    std::sort(lst.begin(), lst.end());

    typedef std::vector<std::pair<int, int> >::iterator vp_iterator;
    std::pair<vp_iterator, vp_iterator> r0;
    std::pair<vp_iterator, vp_iterator> r1;

    for (; first != last; ++first) {
	rp = kmer_fraction(*first);
	if (rp.count > jmin) continue;

	r0 = std::equal_range(lst.begin(), lst.end(), std::make_pair(first->id0, 0), less1st<int, int>);
	r1 = std::equal_range(r0.second, lst.end(), std::make_pair(first->id1, 0), less1st<int, int>);

	first->count += jaz::intersection_size(r0.first, r0.second, r1.first, r1.second, less2nd<int, int>);
    }
} // update_counts2


template <typename Iter>
inline void update_counts(Iter first, Iter last, const std::vector<int>& rem_list, double jmin) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("update_counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    std::vector<std::vector<int>::const_iterator> rem_index;
    jaz::find_all(rem_list.begin(), rem_list.end(), std::back_inserter(rem_index), std::bind2nd(std::equal_to<int>(), -1));

    unsigned int l = rem_index.size();
    read_pair rp;

    for (; first != last; ++first) {
	std::vector<int>::const_iterator begin = rem_list.begin();

	rp = kmer_fraction(*first);
	if (rp.count > jmin) continue;

	for (unsigned int i = 0; i < l; ++i) {
	    if ((first->count + l - i) < (0.01 * jmin * first->size)) break;

	    std::vector<int>::const_iterator end = rem_index[i];

	    if (std::binary_search(begin, end, first->id0) && std::binary_search(begin, end, first->id1)) {
		first->count++;
	    }

	    begin = end + 1;
	} // for i

    } // for first
} // update_counts


void aggregate_rem_list(MPI_Comm comm, std::vector<int>& rem_list) {
#ifdef WITH_MPE
    mpix::MPE_Log mpe_log("aggregate_rem_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // sort rem_list for fast search
    std::vector<std::vector<int>::iterator> rem_index;
    jaz::find_all(rem_list.begin(), rem_list.end(), std::back_inserter(rem_index), std::bind2nd(std::equal_to<int>(), -1));

    std::vector<int>::iterator it = rem_list.begin();

    for (unsigned int i = 0; i < rem_index.size(); ++i) {
	std::sort(it, rem_index[i]);
	it = rem_index[i] + 1;
    }

    // aggregate rem_list
    int rl_sz = rem_list.size();

    std::vector<int> rem_list_sz(size, 0);
    std::vector<int> rem_list_off(size, 0);

    MPI_Allgather(&rl_sz, 1, MPI_INT, &rem_list_sz[0], 1, MPI_INT, comm);

    unsigned int g_rl_sz = std::accumulate(rem_list_sz.begin(), rem_list_sz.end(), 0);
    std::vector<int> rem_list_global(g_rl_sz);

    for (unsigned int i = 1; i < size; ++i) rem_list_off[i] = rem_list_off[i - 1] + rem_list_sz[i - 1];

    MPI_Allgatherv(&rem_list[0], rl_sz, MPI_INT, &rem_list_global[0], &rem_list_sz[0], &rem_list_off[0], MPI_INT, comm);
    rem_list = rem_list_global;
}; // aggregate_rem_list


template <typename Hash>
inline std::pair<bool, std::string> extract_seq_pairs(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
						      std::vector<sketch_id>& sketch_list, Hash hash, std::vector<read_pair>& edges) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // create rem_list and remove singletons
    std::sort(sketch_list.begin(), sketch_list.end());
    std::vector<int> rem_list;

    typedef std::vector<sketch_id>::iterator si_iterator;

    si_iterator iter = sketch_list.begin();
    si_iterator end = sketch_list.end();

    unsigned int max_part = 0;

    while (iter != end) {
	si_iterator temp = jaz::range(iter, end);
	if (temp - iter >= opt.cmax) {
	    // we store in the temporal list (-1 separates lists)
	    for (si_iterator i = iter; i != temp; ++i) {
		if (rem_list.empty() || (rem_list.back() != i->id)) {
		    rem_list.push_back(i->id);
		}
	    }

	    rem_list.push_back(-1);

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

    // balance sketches
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
	unsigned int gmax_part = 0;
	MPI_Allreduce(&max_part, &gmax_part, 1, MPI_UNSIGNED, MPI_MAX, comm);
	part = gmax_part / partition_level(size);
    }

    unsigned int part_id = rank * (std::numeric_limits<unsigned int>::max() / size) + 1;

    sketch_list = mpix::partition_balance(sketch_list, MPI_SKETCH_ID, comm);
    part_id = decompose_sketch_list(sketch_list, part, part_id);

    sketch_list = mpix::partition_balance(sketch_list, MPI_SKETCH_ID, comm);
    part_id = decompose_sketch_list(sketch_list, std::max(part / 4, 500), part_id);

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

#ifdef WITH_MPE
    mpe_log.init("balance sketch_list", "red");
    mpe_log.start();
#endif // WITH_MPE

    sketch_list = mpix::partition_balance(sketch_list, std::equal_to<sketch_id>(),
					  sketch_part_cost<si_iterator>, MPI_SKETCH_ID, 0, comm);

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    MPI_Type_free(&MPI_SKETCH_ID);

#ifdef WITH_MPE
    mpe_log.init("get counts", "red");
    mpe_log.start();
#endif // WITH_MPE

    // get local counts
    std::vector<read_pair> counts;

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

    // perform counts compaction
    std::sort(counts.begin(), counts.end());
    counts.erase(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()), counts.end());

    // perform global reduction
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    counts = mpix::simple_partition(counts, hash, MPI_READ_PAIR, comm);

    std::sort(counts.begin(), counts.end());
    counts.erase(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()), counts.end());

    MPI_Type_free(&MPI_READ_PAIR);

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    // aggregate rem_list
    aggregate_rem_list(comm, rem_list);

    // update counts
    update_counts2(counts.begin(), counts.end(), rem_list, opt.jmin);

#ifdef WITH_MPE
    mpe_log.init("filter candidate edges", "red");
    mpe_log.start();
#endif // WITH_MPE

    // approximate kmer fraction and filter edges
    // count is now approximated kmer fraction
    std::transform(counts.begin(), counts.end(), counts.begin(), kmer_fraction);
    counts.erase(std::remove_if(counts.begin(), counts.end(), not_similar(opt.jmin)), counts.end());

    // store filtered edges (replace with inplace_merge?)
    std::copy(counts.begin(), counts.end(), std::back_inserter(edges));
    std::sort(edges.begin(), edges.end());

    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

#ifdef WITH_MPE
    mpe_log.stop();
#endif // WITH_MPE

    return std::make_pair(true, "");
} // extract_seq_pairs


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

    report << info << "running " << end << " iteration(s): ";

    // main loop (we look for several evidences to get edge)
    for (unsigned int i = 0; i < end; ++i) {
	report << "." << std::flush;

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

#ifdef WITH_MPE
	mpe_log.stop();
#endif // WITH_MPE

#ifdef WITH_MPE
	mpe_log.init("group sketch_list", "red");
	mpe_log.start();
#endif // WITH_MPE

	// group globally sketch list
	sketch_list = mpix::simple_partition(sketch_list, hash_sketch_id, MPI_SKETCH_ID, comm);

#ifdef WITH_MPE
	mpe_log.stop();
#endif // WITH_MPE

	// if (sketch_list.empty()) {
	//     report.critical << warning << "{" << rank << "}" << " empty list of sketches!" << std::endl;
	// }

	// add to edges read pairs with common sketches
	extract_seq_pairs(opt, log, report, comm, sketch_list, hash_read_pair0, edges);
	sketch_list.clear();
    } // for i

    MPI_Type_free(&MPI_SKETCH_ID);

    // rebalance graph for the next stage
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    edges = mpix::simple_partition(edges, hash_read_pair2(SL.N, size), MPI_READ_PAIR, comm);

    MPI_Type_free(&MPI_READ_PAIR);

    report << "" << std::endl;

    return std::make_pair(true, "");
} // generate_edges

#endif // GENERATE_EDGES_HPP
