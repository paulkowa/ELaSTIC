/***
 *  $Id$
 **
 *  File: partition_balance.hpp
 *  Created: Jan 25, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix.
 */

#ifndef MPIX_PARTITION_BALANCE_HPP
#define MPIX_PARTITION_BALANCE_HPP

#include <algorithm>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <mpi.h>


namespace mpix {

  namespace detail {

    struct part {
	int rank;
	unsigned int pos;
	unsigned int S;
    }; // struct part

    inline bool operator<(const part& p1, const part& p2) {
	return (p1.S < p2.S);
    } // operator<

    inline bool block_compare(const part& p1, const part& p2) {
	return (p1.rank < p2.rank) || (!(p2.rank < p1.rank) && p1.pos < p2.pos);
    } // block_compare

    struct local_part {
	explicit local_part(int rank) : rank_(rank) { }
	bool operator()(const part& p) const { return p.rank == rank_; }
	int rank_;
    }; // struct local_part

    inline std::ostream& operator<<(std::ostream& os, const part& p) {
	os << "(" << p.rank << "," << p.pos << "," << p.S << ")";
	return os;
    } // operator<<

    inline void print(const std::vector<std::vector<part> >& list, std::ostream& os) {
	unsigned int k = list.size();
	for (unsigned int i = 0; i < k; ++i) {
	    os << i << ": ";
	    std::copy(list[i].begin(), list[i].end(), std::ostream_iterator<part>(os, " "));
	    os << "\n";
	}
    } // print


    template <typename Fun>
    void balance(std::vector<std::vector<part> >& sched, Fun fun) {
	unsigned int size = sched.size();
	std::vector<unsigned int> S(size, 0);

	std::vector<part> parts;

	for (unsigned int i = 0; i < size; ++i) {
	    std::copy(sched[i].begin(), sched[i].end(), std::back_inserter(parts));
	    sched[i].clear();
	}

	std::sort(parts.begin(), parts.end());
	unsigned int n = parts.size();

	for (unsigned int i = 0; i < n; ++i) {
	    unsigned int pos = std::min_element(S.begin(), S.end()) - S.begin();
	    sched[pos].push_back(parts.back());
	    S[pos] += fun(parts.back().S);
	    parts.pop_back();
	}

    } // balance

    template <typename T> T linear(T x) { return x; }

  } // detail


  template <typename Iter, typename Pred, typename Fun>
  std::pair<typename std::iterator_traits<Iter>::pointer,
	    typename std::iterator_traits<Iter>::pointer>
  partition_balance(Iter first, Iter last, Pred pred, Fun fun, MPI_Datatype Type, int root, MPI_Comm Comm) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      typedef typename std::iterator_traits<Iter>::pointer pointer_type;

      unsigned int n = std::distance(first, last);

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      std::vector<value_type> data(n);
      std::copy(first, last, data.begin());

      // step 1: local sort
      std::sort(data.begin(), data.end(), pred);
      if (size == 1) return std::make_pair<pointer_type, pointer_type>(0, 0);

      // step 2: gather bucket data
      std::vector<unsigned int> parts;

      std::vector<int> parts_sz(size, 0);
      std::vector<int> parts_disp(size, 0);

      std::vector<unsigned int> parts_all;

      unsigned int pos = 0;

      for (unsigned int i = 1; i < n; ++i) {
	  if (pred(data[pos], data[i]) == true) {
	      parts.push_back(i - pos);
	      pos = i;
	  }
      } // for i

      if (n > 0) parts.push_back(n - pos);

      int s = parts.size();
      MPI_Gather(&s, 1, MPI_INT, &parts_sz[0], 1, MPI_INT, root, Comm);

      if (rank == root) {
	  std::partial_sum(parts_sz.begin(), parts_sz.end() - 1, parts_disp.begin() + 1);
	  parts_all.resize(std::accumulate(parts_sz.begin(), parts_sz.end(), 0));
      }

      MPI_Gatherv(&parts[0], s, MPI_UNSIGNED, &parts_all[0], &parts_sz[0], &parts_disp[0], MPI_UNSIGNED, root, Comm);

      // step 3: generate balancing moves
      std::vector<std::vector<detail::part> > moves;

      std::vector<detail::part> moves_buf;
      std::vector<int> moves_sz(size, 0);
      std::vector<int> moves_disp(size, 0);

      if (rank == root) {
	  std::vector<std::vector<detail::part> > sched(size);
	  pos = 0;

	  for (unsigned int i = 0; i < size; ++i) {
	      for (unsigned int j = 0; j < parts_sz[i]; ++j) {
		  detail::part p;
		  p.rank = i;
		  p.pos = j;
		  p.S = parts_all[pos];
		  sched[i].push_back(p);
		  pos++;
	      } // for j
	  } // for i

	  std::vector<unsigned int>().swap(parts_all);

	  detail::balance(sched, fun);
	  moves.resize(size);

	  for (unsigned int i = 0; i < size; ++i) {
	      unsigned int l = sched[i].size();
	      for (unsigned int j = 0; j < l; ++j) {
		  moves[sched[i][j].rank].push_back(sched[i][j]);
		  moves[sched[i][j].rank].back().rank = i;
	      }
	  }
      } // if rank

      // step 4: compact moves
      if (rank == root) {
	  for (unsigned int i = 0; i < size; ++i) {
	      // clean
	      moves[i].erase(std::remove_if(moves[i].begin(), moves[i].end(), detail::local_part(i)), moves[i].end());
	      std::sort(moves[i].begin(), moves[i].end(), detail::block_compare);

	      // merge
	      unsigned int pos = 0;
	      unsigned int res = 0;

	      unsigned int l = moves[i].size();
	      if (l == 0) continue;

	      for (unsigned int j = 1; j < l; ++j) {
		  if ((moves[i][j - 1].rank != moves[i][j].rank) ||
		      (moves[i][j - 1].pos != moves[i][j].pos - 1)) {
		      moves[i][res] = moves[i][pos];
		      for (unsigned int k = pos + 1; k < j; ++k) moves[i][res].S += moves[i][k].S;
		      pos = j;
		      res++;
		  } // if
	      } // for j

	      moves[i][res] = moves[i][pos];
	      for (unsigned int k = pos + 1; k < l; ++k) moves[i][res].S += moves[i][k].S;
	      res++;

	      // clean
	      moves[i].erase(moves[i].begin() + res, moves[i].end());

	      // put into buffer
	      std::copy(moves[i].begin(), moves[i].end(), std::back_inserter(moves_buf));
	      moves_sz[i] = moves[i].size();
	  } // for i

	  std::partial_sum(moves_sz.begin(), moves_sz.end() - 1, moves_disp.begin() + 1);
      } // if rank

      std::vector<std::vector<detail::part> >().swap(moves);

      // Step 5: distribute moves
      int msz = 0;
      MPI_Scatter(&moves_sz[0], 1, MPI_INT, &msz, 1, MPI_INT, root, Comm);

      std::vector<detail::part> my_moves(msz);

      MPI_Datatype MPI_PART;
      MPI_Type_contiguous(sizeof(detail::part), MPI_BYTE, &MPI_PART);
      MPI_Type_commit(&MPI_PART);

      MPI_Scatterv(&moves_buf[0], &moves_sz[0], &moves_disp[0], MPI_PART, &my_moves[0], msz, MPI_PART, root, Comm);

      MPI_Type_free(&MPI_PART);

      // Step 6: prepare data to exchange
      unsigned int l = my_moves.size();
      unsigned int S = 0;

      for (unsigned int i = 0; i < l; ++i) S += my_moves[i].S;

      parts_disp.resize(parts.size());

      if (parts.empty() == false) {
	  parts_disp[0] = 0;
	  std::partial_sum(parts.begin(), parts.end() - 1, parts_disp.begin() + 1);
      }

      std::vector<value_type> send_buf(S);
      std::vector<int> send_sz(size, 0);
      std::vector<int> send_disp(size, 0);

      std::vector<int> send_sz_all(size, 0);
      std::vector<int> send_disp_all(size, 0);

      pos = 0;

      for (unsigned int i = 0; i < l; ++i) {
	  std::copy(data.begin() + parts_disp[my_moves[i].pos],
		    data.begin() + parts_disp[my_moves[i].pos] + my_moves[i].S,
		    send_buf.begin() + pos);
	  send_sz[my_moves[i].rank] += my_moves[i].S;
	  pos += my_moves[i].S;
      }

      std::partial_sum(send_sz.begin(), send_sz.end() - 1, send_disp.begin() + 1);

      for (unsigned int i = 0; i < l; ++i) {
	  pos = parts_disp[my_moves[i].pos];
	  data.erase(data.begin() + parts_disp[my_moves[i].pos],
		     data.begin() + parts_disp[my_moves[i].pos] + my_moves[i].S);
	  for (unsigned int j = 0; j < parts_disp.size(); ++j) {
	      if (pos < parts_disp[j]) parts_disp[j] -= my_moves[i].S;
	  }
      } // for i

      // Step 7: exchange data
      MPI_Alltoall(&send_sz[0], 1, MPI_INT, &send_sz_all[0], 1, MPI_INT, Comm);

      n = data.size();
      S = std::accumulate(send_sz_all.begin(), send_sz_all.end(), 0);

      std::partial_sum(send_sz_all.begin(), send_sz_all.end() - 1, send_disp_all.begin() + 1);

      pointer_type new_data = new value_type[n + S];
      std::copy(data.begin(), data.end(), new_data);

      MPI_Alltoallv(&send_buf[0], &send_sz[0], &send_disp[0], Type,
		    new_data + n, &send_sz_all[0], &send_disp_all[0], Type, Comm);

      return std::make_pair(new_data, new_data + n + S);
  } // partition_balance


  template <typename Iter>
  std::pair<typename std::iterator_traits<Iter>::pointer,
	    typename std::iterator_traits<Iter>::pointer>
  partition_balance(Iter first, Iter last, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return partition_balance(first, last, std::less<value_type>(), detail::linear<unsigned int>, Type, 0, Comm);
  } // partition_balance

} // namespace mpix

#endif // MPIX_PARTITION_BALANCE_HPP
