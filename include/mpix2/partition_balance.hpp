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
 *  This file is part of mpix2.
 */

#ifndef MPIX2_PARTITION_BALANCE_HPP
#define MPIX2_PARTITION_BALANCE_HPP

#include <functional>
#include <vector>
#include <mpi.h>

#include <boost/heap/fibonacci_heap.hpp>


namespace mpix {

  namespace detail {

    struct part {
	int rank;
	int pos;
	int S;
	int cost;
    }; // struct part

    inline bool operator<(const part& p1, const part& p2) {
	return (p1.cost < p2.cost);
    } // operator<

    inline bool block_compare(const part& p1, const part& p2) {
	return (p1.rank < p2.rank) || (!(p2.rank < p1.rank) && (p1.pos < p2.pos));
    } // block_compare

    struct local_part {
	explicit local_part(int rank) : rank_(rank) { }
	bool operator()(const part& p) const { return p.rank == rank_; }
	int rank_;
    }; // struct local_part

    inline std::ostream& operator<<(std::ostream& os, const part& p) {
	os << "(" << p.rank << "," << p.pos << "," << p.S << "," << p.cost << ")";
	return os;
    } // operator<<

    struct balance_heap_node {
	explicit balance_heap_node(int apos = 0, int acost = 0) : pos(apos), cost(acost) { }
	bool operator<(const balance_heap_node& node) const { return node.cost < cost; }
	int pos;
	int cost;
    }; // balance_heap_node

    inline void balance(std::vector<part>& parts, std::vector<std::vector<part> >& sched) {
	int size = sched.size();
	int n = parts.size();

	std::sort(parts.begin(), parts.end());

	typedef boost::heap::fibonacci_heap<balance_heap_node> heap_type;
	std::vector<heap_type::handle_type> Sh(size);
	heap_type S;

	int lim = std::min(size, n);

	// initialize heap
	for (int i = 0; i < lim; ++i) {
	    sched[i].push_back(parts.back());
	    Sh[i] = S.push(balance_heap_node(i, parts.back().cost));
	    parts.pop_back();
	}

	// update heap (and schedule)
	for (int i = 0; i < n - lim; ++i) {
	    balance_heap_node nd = S.top();
	    sched[nd.pos].push_back(parts.back());
	    nd.cost += parts.back().cost;
	    S.update(Sh[nd.pos], nd);
	    parts.pop_back();
	}
    } // balance

    // this function comes from jaz0x
    template <typename Iter, typename Comp>
    Iter range(Iter first, Iter last, Comp comp) {
	if (first == last) return last;
	Iter iter = first;
	for (; (iter != last) && comp(*first, *iter); ++iter);
	return iter;
    } // range

    template <typename Iter> inline int linear(Iter first, Iter last) { return last - first; }

  } // detail


  template <typename Sequence, typename Pred, typename Fun>
  std::vector<typename Sequence::value_type> partition_balance(Sequence& seq,
							       Pred pred,
							       Fun fun,
							       MPI_Datatype Type,
							       int root,
							       MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      typedef typename Sequence::iterator iterator;

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      if (size == 1) {
	  std::vector<value_type> ndata(seq.begin(), seq.end());
	  { Sequence().swap(seq); }
	  return ndata;
      }

      // step 1: gather bucket data
      std::vector<int> parts;
      std::vector<int> parts_all;

      std::vector<int> cost;
      std::vector<int> cost_all;

      std::vector<int> parts_sz(size, 0);
      std::vector<int> parts_disp(size, 0);

      iterator iter = seq.begin();
      iterator end = seq.end();

      while (iter != end) {
	  iterator temp = detail::range(iter, end, pred);
	  parts.push_back(temp - iter);
	  cost.push_back(fun(iter, temp));
	  iter = temp;
      }

      int psz = parts.size();

      MPI_Gather(&psz, 1, MPI_INT, &parts_sz[0], 1, MPI_INT, root, Comm);

      if (rank == root) {
	  std::partial_sum(parts_sz.begin(), parts_sz.end() - 1, parts_disp.begin() + 1);
	  parts_all.resize(std::accumulate(parts_sz.begin(), parts_sz.end(), 0));
      }

      MPI_Gatherv(&parts[0], psz, MPI_INT, &parts_all[0], &parts_sz[0], &parts_disp[0], MPI_INT, root, Comm);

      cost_all.resize(parts_all.size());

      MPI_Gatherv(&cost[0], psz, MPI_INT, &cost_all[0], &parts_sz[0], &parts_disp[0], MPI_INT, root, Comm);

      // step 2: generate balancing moves
      std::vector<int> moves_sz(size, 0);
      std::vector<int> moves_disp(size, 0);

      std::vector<std::vector<detail::part> > moves;

      if (rank == root) {
	  std::vector<std::vector<detail::part> > sched(size);
	  std::vector<detail::part> parts_desc(parts_all.size());

	  int pos = 0;

	  for (int i = 0; i < size; ++i) {
	      for (int j = 0; j < parts_sz[i]; ++j) {
		  parts_desc[pos].rank = i;
		  parts_desc[pos].pos = j;
		  parts_desc[pos].S = parts_all[pos];
		  parts_desc[pos].cost = cost_all[pos];
		  pos++;
	      } // for j
	  } // for i

	  { std::vector<int>().swap(parts_all); }
	  { std::vector<int>().swap(cost_all); }

	  detail::balance(parts_desc, sched);

	  { std::vector<detail::part>().swap(parts_desc); }

	  moves.resize(size);

	  for (int i = 0; i < size; ++i) {
	      int l = sched[i].size();
	      for (int j = 0; j < l; ++j) {
		  moves[sched[i][j].rank].push_back(sched[i][j]);
		  moves[sched[i][j].rank].back().rank = i;
	      }
	  }
      } // if rank

      // step 3: compact moves
      std::vector<detail::part> moves_buf;

      if (rank == root) {
	  for (int i = 0; i < size; ++i) {
	      // clean
	      moves[i].erase(std::remove_if(moves[i].begin(), moves[i].end(), detail::local_part(i)), moves[i].end());
	      std::sort(moves[i].begin(), moves[i].end(), detail::block_compare);

	      // merge
	      int pos = 0;
	      int res = 0;

	      int l = moves[i].size();
	      if (l == 0) continue;

	      for (int j = 1; j < l; ++j) {
		  if ((moves[i][j - 1].rank != moves[i][j].rank) ||
		      (moves[i][j - 1].pos != moves[i][j].pos - 1)) {
		      moves[i][res] = moves[i][pos];
		      for (int k = pos + 1; k < j; ++k) moves[i][res].S += moves[i][k].S;
		      pos = j;
		      res++;
		  } // if
	      } // for j

	      moves[i][res] = moves[i][pos];
	      for (int k = pos + 1; k < l; ++k) moves[i][res].S += moves[i][k].S;
	      res++;

	      // clean
	      moves[i].erase(moves[i].begin() + res, moves[i].end());

	      // put into buffer
	      std::copy(moves[i].begin(), moves[i].end(), std::back_inserter(moves_buf));
	      moves_sz[i] = moves[i].size();
	  } // for i

	  std::partial_sum(moves_sz.begin(), moves_sz.end() - 1, moves_disp.begin() + 1);
      } // if rank

      { std::vector<std::vector<detail::part> >().swap(moves); }

      // Step 4: distribute moves
      int msz = 0;
      MPI_Scatter(&moves_sz[0], 1, MPI_INT, &msz, 1, MPI_INT, root, Comm);

      std::vector<detail::part> my_moves(msz);

      MPI_Datatype MPI_PART;
      MPI_Type_contiguous(sizeof(detail::part), MPI_BYTE, &MPI_PART);
      MPI_Type_commit(&MPI_PART);

      MPI_Scatterv(&moves_buf[0], &moves_sz[0], &moves_disp[0], MPI_PART, &my_moves[0], msz, MPI_PART, root, Comm);

      MPI_Type_free(&MPI_PART);

      // Step 5: prepare data to exchange
      int l = my_moves.size();
      int S = 0;

      for (int i = 0; i < l; ++i) S += my_moves[i].S;

      parts_disp.resize(parts.size());

      if (!parts.empty()) {
	  parts_disp[0] = 0;
	  std::partial_sum(parts.begin(), parts.end() - 1, parts_disp.begin() + 1);
      }

      std::vector<value_type> send_buf(S);

      std::vector<int> send_sz(size, 0);
      std::vector<int> send_disp(size, 0);

      std::vector<int> send_sz_all(size, 0);
      std::vector<int> send_disp_all(size, 0);

      int pos = 0;
      std::vector<bool> erase(seq.size(), false);

      for (int i = 0; i < l; ++i) {
	  std::copy(seq.begin() + parts_disp[my_moves[i].pos],
		    seq.begin() + parts_disp[my_moves[i].pos] + my_moves[i].S,
		    send_buf.begin() + pos);
	  for (int j = parts_disp[my_moves[i].pos]; j < parts_disp[my_moves[i].pos] + my_moves[i].S; ++j) {
	      erase[j] = true;
	  }
	  send_sz[my_moves[i].rank] += my_moves[i].S;
	  pos += my_moves[i].S;
      }

      std::partial_sum(send_sz.begin(), send_sz.end() - 1, send_disp.begin() + 1);

      std::vector<value_type> ndata(seq.size() - pos);
      pos = 0;

      for (int i = 0; i < seq.size(); ++i) if (!erase[i]) ndata[pos++] = seq[i];

      { Sequence().swap(seq); }

      // Step 6: exchange data
      MPI_Alltoall(&send_sz[0], 1, MPI_INT, &send_sz_all[0], 1, MPI_INT, Comm);

      S = std::accumulate(send_sz_all.begin(), send_sz_all.end(), 0);
      std::partial_sum(send_sz_all.begin(), send_sz_all.end() - 1, send_disp_all.begin() + 1);

      int n = ndata.size();
      ndata.resize(n + S);

      MPI_Alltoallv(&send_buf[0], &send_sz[0], &send_disp[0], Type,
		    &ndata[n], &send_sz_all[0], &send_disp_all[0], Type, Comm);

      return ndata;
  } // partition_balance

  template <typename Sequence>
  std::vector<typename Sequence::value_type>
  partition_balance(Sequence& seq, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      typedef typename Sequence::iterator iterator;
      return partition_balance(seq, std::equal_to<value_type>(), detail::linear<iterator>, Type, 0, Comm);
  } // partition_balance

} // namespace mpix

#endif // MPIX2_PARTITION_BALANCE_HPP
