/***
 *  $Id$
 **
 *  File: sample_sort.hpp
 *  Created: Feb 21, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2010-2014 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_SAMPLE_SORT_HPP
#define MPIX2_SAMPLE_SORT_HPP

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <mpi.h>


namespace mpix {

  class bad_sample : public virtual std::exception {
      virtual const char* what() const throw() { return "mpix::bad_sample"; }
  }; // class bad_sample

  namespace detail {

    template <typename Sequence, typename Pred>
    bool sampling(const Sequence& seq, Sequence& pivots, Pred pred, int C,
		  MPI_Datatype Type, MPI_Comm Comm) {
	typedef typename Sequence::value_type value_type;
	typedef typename Sequence::const_iterator iterator;

	int size = 0;
	int rank = 0;

	MPI_Comm_size(Comm, &size);
	MPI_Comm_rank(Comm, &rank);

	long long int n = seq.size();
	int p = size;

	// get total size
	long long int N = 0;
	MPI_Allreduce(&n, &N, 1, MPI_LONG_LONG, MPI_SUM, Comm);

	// if N < 2 nothing to sort
	if (N < 0) throw bad_sample();
	if (N < 2) return false;

	// get sample
	int d = std::max(N / (C * p), 2LL);
	int s = n / d;

	// we always sample smallest and largest element
	Sequence sb(s + 1);

	if (s > 0) for (int i = 0; i < s; ++i) sb[i] = seq[i * d];

	if (n > 0) {
	    sb.back() = seq.back();
	    s++;
	}

	std::vector<int> all_s(p, 0);
	std::vector<int> displ_s(p, 0);

	MPI_Allgather(&s, 1, MPI_INT, &all_s[0], 1, MPI_INT, Comm);

	int S = std::accumulate(all_s.begin(), all_s.end(), 0);
	if (S < 2) throw bad_sample();

	std::partial_sum(all_s.begin(), all_s.end() - 1, displ_s.begin() + 1);

	std::vector<value_type> sample(S);
	MPI_Allgatherv(&sb[0], s, Type, &sample[0], &all_s[0], &displ_s[0], Type, Comm);

	std::sort(sample.begin(), sample.end(), pred);

	// get histogram
	std::vector<int> hist(S, 0);
	iterator iter = seq.begin();

	for (int i = 0; i < S - 1; ++i) {
	    iterator temp = std::upper_bound(iter, seq.end(), sample[i], pred);
	    hist[i] = temp - iter;
	    iter = temp;
	}

	hist.back() = seq.end() - iter;

	std::vector<int> ghist(S, 0);
	MPI_Allreduce(&hist[0], &ghist[0], S, MPI_INT, MPI_SUM, Comm);

	// find pivots
	pivots.resize(p - 1);

	int sz = ghist[0];
	int pos = 0;

	d = std::max(N / p, 3LL);

	// we use greedy approach here
	for (int i = 1; (i < S - 1) && (pos < p - 1); ++i) {
	    if (d < (sz + ghist[i])) {
		if ((d - sz) < sz + ghist[i] - d) {
		    pivots[pos] = sample[i - 1];
		} else {
		    pivots[pos] = sample[i];
		    ++i;
		}
		sz = 0;
		++pos;
	    } // if d
	    sz += ghist[i];
	} // for i

	return true;
    } // sampling

  } // namespace detail


  template <typename Sequence, typename Pred>
  void sample_sort(Sequence& seq, Pred pred, int s, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      typedef typename Sequence::iterator iterator;

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      int n = seq.size();
      int p = size;

      // local sort
      std::sort(seq.begin(), seq.end(), pred);
      if (p == 1) return;

      // pivots
      if (s == 0) s = 16;
      std::vector<value_type> pivots;
      if (detail::sampling(seq, pivots, pred, s, Type, Comm) == false) return;

      // local buffers
      std::vector<int> bin_sz(p, 0);
      std::vector<int> bin_disp(p, 0);

      iterator iter = seq.begin();

      for (int i = 0; i < p - 1; ++i) {
	  iterator temp = std::upper_bound(iter, seq.end(), pivots[i], pred);
	  bin_sz[i] = temp - iter;
	  iter = temp;
      }

      bin_sz.back() = seq.end() - iter;

      std::partial_sum(bin_sz.begin(), bin_sz.end() - 1, bin_disp.begin() + 1);

      // exchange buffers
      std::vector<int> all_bin_sz(p, 0);
      std::vector<int> all_bin_disp(p, 0);

      MPI_Alltoall(&bin_sz[0], 1, MPI_INT, &all_bin_sz[0], 1, MPI_INT, Comm);

      int S = std::accumulate(all_bin_sz.begin(), all_bin_sz.end(), 0);
      std::partial_sum(all_bin_sz.begin(), all_bin_sz.end() - 1, all_bin_disp.begin() + 1);

      // exchange
      Sequence data(S);

      MPI_Alltoallv(&seq[0], &bin_sz[0], &bin_disp[0], Type,
		    &data[0], &all_bin_sz[0], &all_bin_disp[0], Type, Comm);

      seq = data;

      // sort (probably better than merge)
      std::sort(seq.begin(), seq.end(), pred);
  } // sample_sort

  template <typename Sequence>
  void sample_sort(Sequence& seq, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      return sample_sort(seq, std::less<value_type>(), 16, Type, Comm);
  } // sample_sort

} // namespace mpix

#endif // MPIX_SAMPLE_SORT_HPP
