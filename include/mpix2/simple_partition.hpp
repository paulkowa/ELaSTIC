/***
 *  $Id$
 **
 *  File: simple_partition.hpp
 *  Created: Mar 03, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2010-2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_SIMPLE_PARTITION_HPP
#define MPIX2_SIMPLE_PARTITION_HPP

#include <vector>
#include <mpi.h>


namespace mpix {

  template <typename Sequence, typename Hash>
  std::vector<typename Sequence::value_type> simple_partition(Sequence& seq,
							      Hash hash,
							      MPI_Datatype Type,
							      MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;

      int n = seq.size();

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      int p = size;

      if (p == 1) {
	  std::vector<value_type> ndata(seq.begin(), seq.end());
	  { Sequence().swap(seq); }
	  return ndata;
      }

      // analyze local buffer
      std::vector<int> bin_sz(p, 0);
      for (int i = 0; i < n; ++i) bin_sz[hash(seq[i]) % p]++;

      std::vector<int> displ(p, 0);
      for (int i = 1; i < p; ++i) displ[i] = displ[i - 1] + bin_sz[i - 1];

      // prepare send buffer
      std::vector<int> all_bin_sz(p, 0);
      std::vector<value_type> data(n);

      for (typename Sequence::iterator iter(seq.begin()); iter != seq.end(); ++iter) {
	  unsigned int pos = hash(*iter) % p;
	  data[displ[pos] + all_bin_sz[pos]] = *iter;
	  all_bin_sz[pos]++;
      }

      { Sequence().swap(seq); }

      // exchange data
      MPI_Alltoall(&bin_sz[0], 1, MPI_INT, &all_bin_sz[0], 1, MPI_INT, Comm);

      int S = all_bin_sz[0];
      std::vector<int> all_displ(p, 0);

      for (int i = 1; i < p; ++i) {
	  S += all_bin_sz[i];
	  all_displ[i] = all_displ[i - 1] + all_bin_sz[i - 1];
      }

      std::vector<value_type> ndata(S);

      MPI_Alltoallv(&data[0], &bin_sz[0], &displ[0], Type,
		    &ndata[0], &all_bin_sz[0], &all_displ[0], Type, Comm);

      return ndata;
  } // simple_partition

} // namespace mpix

#endif // MPIX2_SIMPLE_PARTITION_HPP
