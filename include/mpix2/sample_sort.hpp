/***
 *  $Id$
 **
 *  File: sample_sort.hpp
 *  Created: Feb 21, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2010-2014 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_SAMPLE_SORT_HPP
#define MPIX2_SAMPLE_SORT_HPP

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <mpi.h>


namespace mpix {

  namespace detail {

    template <typename Sequence>
    void skewed_sampling(const Sequence& seq, Sequence& out,
			 int s, int p,
			 MPI_Datatype Type, int root, MPI_Comm Comm) {
	// p is always > 1
	if (s == 0) s = (0.5 + std::sqrt(p));
	int d = seq.size() / s;

	Sequence sample;

	if (d > 1) {
	    int pos = d >> 1;
	    sample.resize(s);
	    for (int i = 0; i < s; ++i) sample[i] = seq[i * d + pos];
	} else s = 0;

	int rank = 0;
	MPI_Comm_rank(Comm, &rank);

	std::vector<int> all_s(p);
	std::vector<int> displ_s(p, 0);

	MPI_Gather(&s, 1, MPI_INT, &all_s[0], 1, MPI_INT, root, Comm);

	if (rank == root) {
	    int S = std::accumulate(all_s.begin(), all_s.end(), 0);
	    std::partial_sum(all_s.begin(), all_s.end() - 1, displ_s.begin() + 1);
	    out.resize(S);
	} // if rank == 0

	MPI_Gatherv(&sample[0], s, Type, &out[0], &all_s[0], &displ_s[0], Type, root, Comm);
    } // skewed_sampling

  } // namespace detail


  template <typename Sequence, typename Pred>
  void sample_sort(Sequence& seq, int s, Pred pred, MPI_Datatype Type, int root, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      int n = seq.size();
      int p = size;

      // step 1: local sort
      Sequence all_samples;

      // step 2:
      detail::skewed_sampling(seq, all_samples, s, p, Type, root, Comm);

      std::sort(seq.begin(), seq.end(), pred);
      if (p == 1) return;

      // step 3: find and distribute pivots
      std::vector<value_type> sample(p - 1);

      if (rank == root) {
	  std::sort(all_samples.begin(), all_samples.end(), pred);

	  int ss = all_samples.size();
	  if (ss < p - 1) all_samples.resize(p - 1, all_samples[ss - 1]);

	  int d = all_samples.size() / (p - 1);
	  int pos = d >> 1;

	  for (int i = 0; i < p - 1; ++i) sample[i] = all_samples[i * d + pos];
      } // if rank == root

      { Sequence().swap(all_samples); }

      MPI_Bcast(&sample[0], p - 1, Type, root, Comm);

      // step 4: reallocate data
      // step 4a: find bins borders and distribute bins size
      std::vector<int> bins(p, 0);

      // bins stores displacement for local data (bins starting positions)
      for (int i = 1; i < p; ++i) {
	  bins[i] = std::upper_bound(seq.begin(), seq.end(), sample[i - 1], pred) - seq.begin();
      }

      sample.clear();
      std::vector<int> bins_sz(p, 0);

      // bins_sz stores size of each bin
      for (int i = 1; i < p; ++i) {
	  bins_sz[i - 1] = bins[i] - bins[i - 1];
      }
      bins_sz[p - 1] = n - bins[p - 1];

      // all_bins_sz stores size of bins after data exchange
      std::vector<int> all_bins_sz(p, 0);
      MPI_Alltoall(&bins_sz[0], 1, MPI_INT, &all_bins_sz[0], 1, MPI_INT, Comm);

      // step 4b: allocate memory and move data
      int S = std::accumulate(all_bins_sz.begin(), all_bins_sz.end(), 0);

      std::vector<int> all_bins(p, 0);
      std::partial_sum(all_bins_sz.begin(), all_bins_sz.end() - 1, all_bins.begin() + 1);

      Sequence data(S);

      MPI_Alltoallv(&seq[0], &bins_sz[0], &bins[0], Type,
		    &data[0], &all_bins_sz[0], &all_bins[0], Type, Comm);

      seq = data;

      // step 5: merge received data into final sorted range
      // we use sorting instead of merge to save memory
      // extra cost is negligible or can be actually better
      // if we are short of memory
      std::sort(seq.begin(), seq.end(), pred);
  } // sample_sort

  template <typename Sequence, typename Pred>
  void sample_sort(Sequence& seq, Pred pred, MPI_Datatype Type, MPI_Comm Comm) {
      return sample_sort(seq, 0, pred, Type, 0, Comm);
  } // sample_sort

  template <typename Sequence>
  void sample_sort(Sequence& seq, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      return sample_sort(seq, 0, std::less<value_type>(), Type, 0, Comm);
  } // sample_sort

} // namespace mpix

#endif // MPIX_SAMPLE_SORT_HPP
