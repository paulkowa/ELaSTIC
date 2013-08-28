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

#include <numeric>
#include <set>
#include <string>
#include <vector>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <jaz/algorithm.hpp>
#include <mpix/partition_balance.hpp>
#include <mpix/simple_partition.hpp>

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"

#ifdef WITH_MPE
#include <mpe.h>
#endif // WITH_MPE


template <typename T> inline T sqr(T x) { return x * x; }


template <typename Iter>
inline void update_counts(Iter first, Iter last, const std::vector<int>& rem_list, double jmin) {
#ifdef WITH_MPE
    int mpe_start = MPE_Log_get_event_number();
    int mpe_stop = MPE_Log_get_event_number();
    MPE_Describe_state(mpe_start, mpe_stop, "update_counts", "white");
    MPE_Log_event(mpe_start, 0, "start update_counts");
#endif // WITH_MPE

    std::vector<std::vector<int>::const_iterator> rem_index;
    jaz::find_all(rem_list.begin(), rem_list.end(), std::back_inserter(rem_index), std::bind2nd(std::equal_to<int>(), -1));

    unsigned int l = rem_index.size();

    for (; first != last; ++first) {
	std::vector<int>::const_iterator begin = rem_list.begin();

	for (unsigned int i = 0; i < l; ++i) {
	    if ((first->count + l - i) < (0.01 * jmin * first->size)) break;

	    std::vector<int>::const_iterator end = rem_index[i];

	    if (std::binary_search(begin, end, first->id0) && std::binary_search(begin, end, first->id1)) {
		first->count++;
	    }

	    begin = end + 1;
	} // for i

    } // for first

#ifdef WITH_MPE
    MPE_Log_event(mpe_stop, 0, "stop update_counts");
#endif // WITH_MPE
} // update_counts


void aggregate_rem_list(MPI_Comm comm, std::vector<int>& rem_list) {
#ifdef WITH_MPE
    int mpe_start = MPE_Log_get_event_number();
    int mpe_stop = MPE_Log_get_event_number();
    MPE_Describe_state(mpe_start, mpe_stop, "aggregate_rem_list", "white");
    MPE_Log_event(mpe_start, 0, "start aggregate_rem_list");
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

#ifdef WITH_MPE
    MPE_Log_event(mpe_stop, 0, "stop aggregate_rem_list");
#endif // WITH_MPE
}; // aggregate_rem_list


template <typename Hash>
inline std::pair<bool, std::string> extract_seq_pairs(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
						      sketch_id*& sr, sketch_id*& sr_end, Hash hash,
						      std::vector<read_pair>& edges) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // create rem_list
    std::sort(sr, sr_end);
    std::vector<int> rem_list;

    unsigned int pos = 0;
    unsigned int n = (sr_end - sr);

    for (unsigned int i = 1; i < n; ++i) {
	if ((sr[pos] != sr[i]) || (i == n - 1)) {
	    unsigned int end = (i == (n - 1)) ? n : i;

	    if ((end - pos) >= opt.cmax) {
		// we store in the temporal list (-1 separates lists)
		for (unsigned int j = pos; j < end; ++j) {
		    if ((rem_list.empty() == true) || (rem_list.back() != sr[j].id)) {
			rem_list.push_back(sr[j].id);
		    }
		}
		rem_list.push_back(-1);

		// and we remove from sr
		std::copy(sr + i, sr + n, sr + pos);
		n = (i == (n - 1)) ? pos + n - end : pos + n - i;
		i = pos;
	    } else pos = i;
	} // if
    } // for i

    sr_end = sr + n;

    // balance sketches
    MPI_Datatype MPI_SKETCH_ID;
    MPI_Type_contiguous(sizeof(sketch_id), MPI_BYTE, &MPI_SKETCH_ID);
    MPI_Type_commit(&MPI_SKETCH_ID);

    sketch_id* sr_temp = 0;
    boost::tie(sr_temp, sr_end) = mpix::partition_balance(sr, sr_end, sketch_compare, sqr<unsigned int>, MPI_SKETCH_ID, 0, comm);

    MPI_Type_free(&MPI_SKETCH_ID);

    delete[] sr;
    sr = sr_temp;

    std::sort(sr, sr_end);

#ifdef WITH_MPE
    int mpe_gc_start = MPE_Log_get_event_number();
    int mpe_gc_stop = MPE_Log_get_event_number();
    MPE_Describe_state(mpe_gc_start, mpe_gc_stop, "get_counts", "white");
    MPE_Log_event(mpe_gc_start, 0, "start get counts");
#endif // WITH_MPE

    // get local counts
    std::vector<read_pair> counts;
    sr_temp = sr;

    while (sr_temp != sr_end) {
	sketch_id* sr_iter = jaz::range(sr_temp, sr_end); // ???

	// enumerate all pairs with the same sketch
	for (sketch_id* j = sr_temp; j != sr_iter - 1; ++j) {
	    for (sketch_id* k = j + 1; k != sr_iter; ++k) {
		if (j->id != k->id) counts.push_back(make_read_pair(*j, *k));
	    } // for k
	} // for j

	sr_temp = sr_iter;
    } // while

    // perform counts compaction
    std::sort(counts.begin(), counts.end());
    counts.erase(jaz::compact(counts.begin(), counts.end(), std::plus<read_pair>()), counts.end());

    // perform global reduction
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    read_pair* first = 0;
    read_pair* last = 0;

    boost::tie(first, last) = mpix::simple_partition(counts.begin(), counts.end(), hash, MPI_READ_PAIR, 0, comm);
    std::vector<read_pair>().swap(counts);

    std::sort(first, last);
    last = jaz::compact(first, last, std::plus<read_pair>());

    MPI_Type_free(&MPI_READ_PAIR);

#ifdef WITH_MPE
    MPE_Log_event(mpe_gc_stop, 0, "stop get counts");
#endif // WITH_MPE

    // aggregate rem_list
    aggregate_rem_list(comm, rem_list);

    // update counts
    update_counts(first, last, rem_list, opt.jmin);

#ifdef WITH_MPE
    int mpe_fce_start = MPE_Log_get_event_number();
    int mpe_fce_stop = MPE_Log_get_event_number();
    MPE_Describe_state(mpe_fce_start, mpe_fce_stop, "filter candidate edges", "white");
    MPE_Log_event(mpe_fce_start, 0, "start filter candidate edges");
#endif // WITH_MPE

    // approximate kmer fraction and filter edges
    // count is now approximated kmer fraction
    std::transform(first, last, first, kmer_fraction);
    last = std::remove_if(first, last, not_similar(opt.jmin));

    // store filtered edges (replace with inplace_merge?)
    std::copy(first, last, std::back_inserter(edges));
    std::sort(edges.begin(), edges.end());

    delete[] first;

    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

#ifdef WITH_MPE
    MPE_Log_event(mpe_fce_stop, 0, "stop filter candidate edges");
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
	int mpe_start = MPE_Log_get_event_number();
	int mpe_stop = MPE_Log_get_event_number();
	MPE_Describe_state(mpe_start, mpe_stop, "extract sketches", "white");
	MPE_Log_event(mpe_start, 0, "start extract sketches pairs");
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
	MPE_Log_event(mpe_stop, 0, "stop extract sketches");
#endif // WITH_MPE

	// sort globally sketch list (SR)
	sketch_id* first = 0;
	sketch_id* last = 0;

	boost::tie(first, last) = mpix::simple_partition(sketch_list.begin(), sketch_list.end(), hash_sketch_id,
							 MPI_SKETCH_ID, 0, comm);

	sketch_list.clear();

	if (first == last) {
	    report.critical << warning << "{" << rank << "}" << " empty list of sketches!" << std::endl;
	}

	// add to edges read pairs with common sketches
	extract_seq_pairs(opt, log, report, comm, first, last, hash_read_pair2(SL.N, size), edges);

	delete[] first;
    } // for i

    MPI_Type_free(&MPI_SKETCH_ID);
    report << "" << std::endl;

    return std::make_pair(true, "");
} // generate_edges

#endif // GENERATE_EDGES_HPP
