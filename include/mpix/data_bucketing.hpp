/***
 *  $Id$
 **
 *  File: data_bucketing.hpp
 *  Created: Mar 03, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2010 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix.
 */

#ifndef DATA_BUCKETING_HPP
#define DATA_BUCKETING_HPP

#include <iterator>
#include <utility>
#include <vector>
#include <mpi.h>


namespace mpix {

  template <typename Iter, typename Hash>
  std::pair<typename std::iterator_traits<Iter>::pointer,
	    typename std::iterator_traits<Iter>::pointer>
  data_bucketing(Iter first, Iter last, Hash hash, MPI_Datatype Type, int root, MPI_Comm Comm) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;

      unsigned int n = std::distance(first, last);

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      unsigned int p = size;

      if (p == 1) {
	  value_type* data = new value_type[n];
	  std::copy(first, last, data);
	  return std::make_pair(data, data + n);
      }

      std::vector<int> bin_sz(p, 0);
      for (Iter iter = first; iter != last; ++iter) bin_sz[hash(*iter) % p]++;

      std::vector<int> displ(p, 0);
      for (unsigned int i = 1; i < p; ++i) displ[i] = displ[i - 1] + bin_sz[i - 1];

      std::vector<int> all_bin_sz(p, 0);
      value_type* data = new value_type[n];

      for (Iter iter = first; iter != last; ++iter) {
	  unsigned int pos = hash(*iter) % p;
	  data[displ[pos] + all_bin_sz[pos]] = *iter;
	  all_bin_sz[pos]++;
      }

      MPI_Alltoall(&bin_sz[0], 1, MPI_INT, &all_bin_sz[0], 1, MPI_INT, Comm);

      unsigned int S = all_bin_sz[0];
      std::vector<int> all_displ(p, 0);

      for (unsigned int i = 1; i < p; ++i) {
	  S += all_bin_sz[i];
	  all_displ[i] = all_displ[i - 1] + all_bin_sz[i - 1];
      }

      value_type* new_data = new value_type[S];

      MPI_Alltoallv(data, &bin_sz[0], &displ[0], Type,
		    new_data, &all_bin_sz[0], &all_displ[0], Type, Comm);

      delete[] data;

      return std::make_pair(new_data, new_data + S);
  } // data_bucketing

} // namespace mpix

#endif // DATA_BUCKETING_HPP
