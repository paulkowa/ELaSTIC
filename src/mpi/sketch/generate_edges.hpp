/***
 *  $Id$
 **
 *  File: generate_edges.hpp
 *  Created: May 24, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
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

#include <mpix/data_bucketing.hpp>

#include "SequenceDB.hpp"
#include "config_log.hpp"


template <typename Iter, typename Oper>
Iter compact(Iter first, Iter last, Oper op) {
    if (first == last) return last;

    Iter res = first;
    Iter pos = first;

    first++;

    for (; first != last; ++first) {
	if (*pos < *first) {
	    *res = *pos;
	    pos++;
	    *res = std::accumulate(pos, first, *res, op);
	    pos = first;
	    res++;
	}
    }

    *res = *pos;
    pos++;
    *res = std::accumulate(pos, last, *res, op);
    res++;

    return res;
} // compact


template <typename Iter>
inline void update_counts(Iter first, Iter last,
			  const std::vector<unsigned int>& rem_list,
			  const std::vector<unsigned int>& rem_index) {
    for (; first != last; ++first) {
	// 0 marks end of index list
	unsigned int pos = (std::find(rem_index.begin() + 1, rem_index.end(), 0) - rem_index.begin());

	for (unsigned int i = 1; i < pos; ++i) {
	    std::vector<unsigned int>::const_iterator begin = rem_list.begin() + rem_index[i - 1];
	    std::vector<unsigned int>::const_iterator end = rem_list.begin() + rem_index[i];
	    if (std::binary_search(begin, end, first->id0) && std::binary_search(begin, end, first->id1)) {
		first->count++;
	    }
	} // for i

    } // for first
} // update_counts


inline std::pair<bool, std::string> extract_seq_pairs(const AppConfig& opt, AppLog& log, Reporter& report,
						      const sketch_id* sr, const sketch_id* sr_end,
						      MPI_Comm comm,
						      std::vector<read_pair>& edges) {
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    unsigned int n = (sr_end - sr);

    std::vector<read_pair> counts;

    std::vector<unsigned int> rem_list;
    std::vector<unsigned int> rem_index(1, 0);

    unsigned int pos = 0;

    // get local counts
    for (unsigned int i = 1; i < n; ++i) {
	if ((sr[pos] != sr[i]) || (i == n - 1)) {
	    unsigned int end = (sr[pos] == sr[i]) ? n : i;
	    if ((end - pos) < opt.cmax) {
		// enumerate all pairs with the same sketch
		for (unsigned int j = pos; j < end - 1; ++j) {
		    for (unsigned int k = j + 1; k < end; ++k) {
			counts.push_back(make_read_pair(sr[j], sr[k]));
		    } // for k
		} // for j
	    } else {
		// we store in the temporal list
		for (unsigned int j = pos; j < end; ++j) rem_list.push_back(sr[j].id);
		rem_index.push_back(rem_list.size());
	    }
	    pos = i;
	} // if
    } // for i

    // sort rem_list for fast search
    for (unsigned int i = 1; i < rem_index.size(); ++i) {
	std::sort(rem_list.begin() + rem_index[i - 1], rem_list.begin() + rem_index[i]);
    }

    // perform counts compaction
    std::sort(counts.begin(), counts.end());
    counts.erase(compact(counts.begin(), counts.end(), std::plus<read_pair>()), counts.end());

    // perform global reduction
    MPI_Datatype MPI_READ_PAIR;
    MPI_Type_contiguous(sizeof(read_pair), MPI_BYTE, &MPI_READ_PAIR);
    MPI_Type_commit(&MPI_READ_PAIR);

    read_pair* first = 0;
    read_pair* last = 0;

    // boost::tie(first, last) = mpix::data_bucketing(counts.begin(), counts.end(), hash_read_pair0,
    // 						   MPI_READ_PAIR, 0, comm);

    boost::tie(first, last) = mpix::data_bucketing(counts.begin(), counts.end(), hash_read_pair2(log.input, size),
						   MPI_READ_PAIR, 0, comm);

    std::vector<read_pair>().swap(counts);

    std::sort(first, last);
    last = compact(first, last, std::plus<read_pair>());

    MPI_Type_free(&MPI_READ_PAIR);

    // find size of the largest rem_list
    unsigned int rl_size[2];
    unsigned int rl_size_in[2] = { 0, 0 };

    rl_size[0] = rem_list.size();
    rl_size[1] = rem_index.size();

    MPI_Allreduce(rl_size, rl_size_in, 2, MPI_UNSIGNED, MPI_MAX, comm);

    rem_list.resize(rl_size_in[0], 0);
    rem_index.resize(rl_size_in[1], 0);

    // update candidate pair counts (we need to shift list of updates)
    const int MPI_SHIFT_TAG = 111;

    MPI_Status stat;

    int to = (rank + 1) % size;
    int from = (size + rank - 1) % size;

    std::vector<unsigned int> rem_list_recv(rl_size_in[0], 0);
    std::vector<unsigned int> rem_index_recv(rl_size_in[1], 0);

    // we shift p times list in a circular fashion
    for (unsigned int i = 0; i < size; ++i) {
	update_counts(first, last, rem_list, rem_index);

	MPI_Sendrecv(&rem_list[0], rem_list.size(), MPI_UNSIGNED, to, MPI_SHIFT_TAG,
		     &rem_list_recv[0], rem_list_recv.size(), MPI_UNSIGNED, from, MPI_SHIFT_TAG,
		     comm, &stat);

	MPI_Sendrecv(&rem_index[0], rem_index.size(), MPI_UNSIGNED, to, MPI_SHIFT_TAG,
	 	     &rem_index_recv[0], rem_index_recv.size(), MPI_UNSIGNED, from, MPI_SHIFT_TAG,
		     comm, &stat);

	rem_list.swap(rem_list_recv);
	rem_index.swap(rem_index_recv);
    } // for i

    // approximate Jaccard index and filter edges
    // count is now approximated JI
    std::transform(first, last, first, approximate_jaccard(opt.kmer, opt.mod));
    last = std::remove_if(first, last, not_similar(opt.jmin));

    // store filtered edges (replace with inplace_merge?)
    std::copy(first, last, std::back_inserter(edges));
    std::sort(edges.begin(), edges.end());

    delete[] first;

    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    return std::make_pair(true, "");
} // extract_seq_pairs


inline std::pair<bool, std::string> generate_edges(const AppConfig& opt, AppLog& log, Reporter& report,
						   const std::vector<Sequence>& seqs,
						   const std::vector<shingle_list_type>& shingles,
						   MPI_Comm comm,
						   std::vector<read_pair>& edges) {
    report << step << "extracting candidate edges..." << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    std::vector<sketch_id> sketch_list;

    unsigned int n = seqs.size();
    unsigned int end = std::min(opt.iter, opt.mod);

    // register MPI type for sketch_id pair
    MPI_Datatype MPI_SKETCH_ID;
    MPI_Type_contiguous(sizeof(sketch_id), MPI_BYTE, &MPI_SKETCH_ID);
    MPI_Type_commit(&MPI_SKETCH_ID);

    report << info << "running " << end << " iterations: ";

    // main loop (we look for several evidences to get edge)
    for (unsigned int i = 0; i < end; ++i) {
	report << "." << std::flush;

	// get sketch-read pairs for given mod value
	for (unsigned int j = 0; j < n; ++j) {
	    unsigned int id = seqs[j].id;
	    unsigned short int len = seqs[j].s.size();

	    unsigned int l = shingles[j].size();

	    for (unsigned int k = 0; k < l; ++k) {
		if (shingles[j][k] % opt.mod == i) {
		    sketch_list.push_back(make_sketch_id(shingles[j][k], id, len));
		}
	    }
	} // for j

	// sort globally sketch list (SR)
	sketch_id* first = 0;
	sketch_id* last = 0;

	boost::tie(first, last) = mpix::data_bucketing(sketch_list.begin(), sketch_list.end(), sketch_id_hash,
						       MPI_SKETCH_ID, 0, comm);

	sketch_list.clear();
	std::sort(first, last);

	if (first == last) {
	    report.critical << warning << "{" << rank << "}" << " empty list of sketches!" << std::endl;
	}

	// extract read pairs with common sketches
	extract_seq_pairs(opt, log, report, first, last, comm, edges);

	delete[] first;
    } // for i

    MPI_Type_free(&MPI_SKETCH_ID);

    report << "" << std::endl;

    // get total number of candidate edges
    unsigned int loc = edges.size();

    unsigned int min = 0;
    unsigned int max = 0;
    unsigned int tot = 0;

    MPI_Reduce(&loc, &tot, 1, MPI_UNSIGNED, MPI_SUM, 0, comm);
    MPI_Reduce(&loc, &min, 1, MPI_UNSIGNED, MPI_MIN, 0, comm);
    MPI_Reduce(&loc, &max, 1, MPI_UNSIGNED, MPI_MAX, 0, comm);

    report << info << "found " << tot << " candidate edges" << std::endl;
    report << info << "edges distribution: [" << min << "," << max << "]" << std::endl;

    log.cedges = tot;

    return std::make_pair(true, "");
} // generate_edges

#endif // GENERATE_EDGES_HPP
