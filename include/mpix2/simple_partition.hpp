/***
 *  $Id$
 **
 *  File: simple_partition.hpp
 *  Created: Mar 03, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2010-2014 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_SIMPLE_PARTITION_HPP
#define MPIX2_SIMPLE_PARTITION_HPP

#include <numeric>
#include <vector>
#include <mpi.h>


namespace mpix {

  template <typename Sequence, typename Hash>
  void simple_partition(Sequence& seq, Hash hash, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      int p = size;
      if (p == 1) return;

      int n = seq.size();

      // analyze local buffer
      std::vector<int> bin_sz(p, 0);
      for (int i = 0; i < n; ++i) bin_sz[hash(seq[i]) % p]++;

      std::vector<int> displ(p, 0);
      std::partial_sum(bin_sz.begin(), bin_sz.end() - 1, displ.begin() + 1);

      // prepare send buffer
      std::vector<int> all_bin_sz(p, 0);
      std::vector<value_type> data(n);

      for (int i = 0; i < n; ++i) {
	  int pos = hash(seq[i]) % p;
	  data[displ[pos] + all_bin_sz[pos]] = seq[i];
	  all_bin_sz[pos]++;
      }

      { Sequence().swap(seq); }

      // exchange data
      MPI_Alltoall(&bin_sz[0], 1, MPI_INT, &all_bin_sz[0], 1, MPI_INT, Comm);

      int S = std::accumulate(all_bin_sz.begin(), all_bin_sz.end(), 0);

      std::vector<int> all_displ(p, 0);
      std::partial_sum(all_bin_sz.begin(), all_bin_sz.end() - 1, all_displ.begin() + 1);

      seq.resize(S);

      MPI_Alltoallv(&data[0], &bin_sz[0], &displ[0], Type,
		    &seq[0], &all_bin_sz[0], &all_displ[0], Type, Comm);
  } // simple_partition

} // namespace mpix

#endif // MPIX2_SIMPLE_PARTITION_HPP
