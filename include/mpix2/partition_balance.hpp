/***
 *  $Id$
 **
 *  File: partition_balance.hpp
 *  Created: Jan 25, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2013-2014 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_PARTITION_BALANCE_HPP
#define MPIX2_PARTITION_BALANCE_HPP

#include <algorithm>
#include <functional>
#include <vector>
#include <mpi.h>

#include <boost/heap/fibonacci_heap.hpp>

namespace mpix {

  namespace detail {

    // this function comes from jaz0x
    template <typename Iter, typename Comp>
    Iter range(Iter first, Iter last, Comp comp) {
	if (first == last) return last;
	Iter iter = first;
	for (; (iter != last) && comp(*first, *iter); ++iter);
	return iter;
    } // range


    struct task {
	int size;
	int cost;
    }; // struct task

    inline bool operator<(const task& t1, const task& t2) {
	return (t1.cost < t2.cost);
    } // operator<

    inline task make_task(int size, int cost) {
	task t;
	t.size = size;
	t.cost = cost;
	return t;
    } // make_task

    inline std::ostream& operator<<(std::ostream& os, const task& t) {
	os << "(" << t.size << "," << t.cost << ")";
	return os;
    } // operator<<


    struct ext_task {
	int rank;
	int pos;
	task t;
    }; // struct

    inline ext_task operator+(const ext_task& t1, const ext_task& t2) {
	ext_task t(t1);
	t.t.size += t2.t.size;
	return t;
    } // operator+

    inline bool operator<(const ext_task& t1, const ext_task& t2) {
	return (t1.t < t2.t);
    } // operator<

    inline bool block_compare(const ext_task& t1, const ext_task& t2) {
	return (t1.rank < t2.rank) || (!(t2.rank < t1.rank) && (t1.pos < t2.pos));
    } // block_compare

    inline bool is_task_block(const ext_task& t1, const ext_task& t2) {
	return (t1.rank != t2.rank) || (t1.pos != (t2.pos - 1));
    } // is_task_block

    struct is_local_task {
	explicit is_local_task(int rank) : rank_(rank) { }
	bool operator()(const ext_task& t) const { return t.rank == rank_; }
	int rank_;
    }; // struct is_local_task

    inline ext_task make_ext_task(int rank, int pos, const task& t) {
	ext_task et;
	et.rank = rank;
	et.pos = pos;
	et.t = t;
	return et;
    } // make_ext_task

    inline std::ostream& operator<<(std::ostream& os, const ext_task& t) {
	os << "(" << t.rank << "," << t.pos << "," << t.t << ")";
	return os;
    } // operator<<


    struct heap_node {
	explicit heap_node(int apos = 0, int acost = 0) : pos(apos), cost(acost) { }
	int pos;
	int cost;
    }; // heap_node

    inline bool operator<(const heap_node& n1, const heap_node& n2) {
	return n2.cost < n1.cost;
    } // operator<


    inline void balance(std::vector<ext_task>& tasks, std::vector<std::vector<ext_task> >& sched) {
	int size = sched.size();
	int n = tasks.size();

	std::sort(tasks.begin(), tasks.end());

	typedef boost::heap::fibonacci_heap<heap_node> heap_type;
	std::vector<heap_type::handle_type> Sh(size);
	heap_type S;

	int lim = std::min(size, n);

	// initialize heap
	for (int i = 0; i < lim; ++i) {
	    sched[i].push_back(tasks.back());
	    Sh[i] = S.push(heap_node(i, tasks.back().t.cost));
	    tasks.pop_back();
	}

	// update heap (and schedule)
	for (int i = 0; i < n - lim; ++i) {
	    heap_node nd = S.top();
	    sched[nd.pos].push_back(tasks.back());
	    nd.cost += tasks.back().t.cost;
	    S.update(Sh[nd.pos], nd);
	    tasks.pop_back();
	}
    } // balance


    template <typename Iter> inline int linear(Iter first, Iter last) { return last - first; }

  } // namespace detail


  template <typename Sequence, typename Pred, typename Fun>
  void partition_balance(Sequence& seq, Pred pred, Fun fun,
			 MPI_Datatype Type, int root, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      typedef typename Sequence::iterator iterator;

      int size = 0;
      int rank = 0;

      MPI_Comm_size(Comm, &size);
      MPI_Comm_rank(Comm, &rank);

      if (size == 1) return;

      // step 1: gather tasks data
      std::vector<detail::task> tasks;
      std::vector<detail::task> tasks_all;

      iterator iter = seq.begin();
      iterator end = seq.end();

      while (iter != end) {
	  iterator temp = detail::range(iter, end, pred);
	  tasks.push_back(detail::make_task(temp - iter, fun(iter, temp)));
	  iter = temp;
      } // while

      // we send tasks as 2 integers
      // hence size must be multiplied by 2
      int tsz = 2 * tasks.size();

      std::vector<int> tasks_sz(size, 0);
      std::vector<int> tasks_disp(size, 0);

      MPI_Gather(&tsz, 1, MPI_INT, &tasks_sz[0], 1, MPI_INT, root, Comm);

      if (rank == root) {
	  std::partial_sum(tasks_sz.begin(), tasks_sz.end() - 1, tasks_disp.begin() + 1);
	  tasks_all.resize(std::accumulate(tasks_sz.begin(), tasks_sz.end(), 0) >> 1);
      }

      MPI_Gatherv(reinterpret_cast<int*>(&tasks[0]), tsz, MPI_INT, reinterpret_cast<int*>(&tasks_all[0]), &tasks_sz[0], &tasks_disp[0], MPI_INT, root, Comm);

      // step 2: generate balancing moves
      std::vector<std::vector<detail::ext_task> > moves(size);

      if (rank == 0) {
	  // we extend task information here to minimize message size
	  std::vector<detail::ext_task> ext_tasks_all(tasks_all.size());
	  int pos = 0;

	  for (int i = 0; i < size; ++i) {
	      for (int j = 0; j < (tasks_sz[i] >> 1); ++j, ++pos) {
		  ext_tasks_all[pos] = make_ext_task(i, j, tasks_all[pos]);
	      }
	  } // for i

	  { std::vector<detail::task>().swap(tasks_all); }

	  std::vector<std::vector<detail::ext_task> > sched(size);
	  detail::balance(ext_tasks_all, sched);

	  { std::vector<detail::ext_task>().swap(ext_tasks_all); }

	  // who sends where
	  for (int i = 0; i < size; ++i) {
	      int l = sched[i].size();
	      for (int j = 0; j < l; ++j) {
		  moves[sched[i][j].rank].push_back(sched[i][j]);
		  moves[sched[i][j].rank].back().rank = i;
	      }
	  } // for i
      } // if rank

      std::vector<int> moves_sz(size, 0);
      std::vector<int> moves_disp(size, 0);

      std::vector<detail::ext_task> moves_buf;

      // step 3: compact moves
      typedef std::vector<detail::ext_task>::iterator task_iterator;

      if (rank == 0) {
	  // clean
	  for (int i = 0; i < size; ++i) {
	      moves[i].erase(std::remove_if(moves[i].begin(), moves[i].end(), detail::is_local_task(i)), moves[i].end());
	  }

	  // merge
	  for (int i = 0; i < size; ++i) {
	      if (moves[i].empty()) continue;

	      std::sort(moves[i].begin(), moves[i].end(), detail::block_compare);

	      task_iterator res = moves[i].begin();
	      task_iterator iter = moves[i].begin();
	      task_iterator end = moves[i].end();

	      while (iter < end) {
		  task_iterator temp = std::adjacent_find(iter, end, detail::is_task_block);
		  *(res++) = std::accumulate(iter + 1, std::min(end, temp + 1), *iter);
		  iter = temp + 1;
	      } // while

	      moves[i].erase(res, moves[i].end());
	  } // for i

	  // prepare send buffer
	  for (int i = 0; i < size; ++i) {
	      std::copy(moves[i].begin(), moves[i].end(), std::back_inserter(moves_buf));
	      moves_sz[i] = moves[i].size();
	  }

	  std::partial_sum(moves_sz.begin(), moves_sz.end() - 1, moves_disp.begin() + 1);
      } // if rank

      { std::vector<std::vector<detail::ext_task> >().swap(moves); }

      // step 4: distribute moves
      int msz = 0;
      MPI_Scatter(&moves_sz[0], 1, MPI_INT, &msz, 1, MPI_INT, root, Comm);

      std::vector<detail::ext_task> my_moves(msz);

      MPI_Datatype MPI_EXT_TASK;
      MPI_Type_contiguous(sizeof(detail::ext_task), MPI_BYTE, &MPI_EXT_TASK);
      MPI_Type_commit(&MPI_EXT_TASK);

      MPI_Scatterv(&moves_buf[0], &moves_sz[0], &moves_disp[0], MPI_EXT_TASK, &my_moves[0], msz, MPI_EXT_TASK, root, Comm);

      MPI_Type_free(&MPI_EXT_TASK);

      // step 5: prepare data to exchange
      int S = my_moves.empty() ? 0 : std::accumulate(my_moves.begin() + 1, my_moves.end(), *my_moves.begin()).t.size;
      std::vector<value_type> send_buf(S);

      // get data localization
      tasks_disp.resize(tasks.size());

      if (!tasks.empty()) {
	  tasks_disp[0] = 0;
	  for (int i = 1; i < tasks.size(); ++i) tasks_disp[i] = tasks_disp[i - 1] + tasks[i - 1].size;
      }

      // place data to send
      std::fill(tasks_sz.begin(), tasks_sz.end(), 0);
      iterator out = send_buf.begin();

      // we use erase markers for performance
      std::vector<bool> erase(seq.size(), false);

      for (task_iterator i = my_moves.begin(); i != my_moves.end(); ++i) {
	  int pos = tasks_disp[i->pos];

	  iter = seq.begin() + pos;
	  end = iter + i->t.size;

	  std::fill(erase.begin() + pos, erase.begin() + pos + i->t.size, true);
	  std::copy(iter, end, out);

	  tasks_sz[i->rank] += i->t.size;
	  out += i->t.size;
      } // for i

      // clean local memory (ugly implementation)
      int pos = std::find(erase.begin(), erase.end(), true) - erase.begin();
      int sub = std::find(erase.begin() + pos, erase.end(), false) - erase.begin();

      for (; sub < seq.size(); ++sub) if (!erase[sub]) seq[pos++] = seq[sub];

      // step 6: exchange data
      std::fill(tasks_disp.begin(), tasks_disp.end(), 0);
      std::partial_sum(tasks_sz.begin(), tasks_sz.end() - 1, tasks_disp.begin() + 1);

      std::vector<int> recv_sz(size, 0);
      std::vector<int> recv_disp(size, 0);

      MPI_Alltoall(&tasks_sz[0], 1, MPI_INT, &recv_sz[0], 1, MPI_INT, Comm);

      S = std::accumulate(recv_sz.begin(), recv_sz.end(), 0);
      std::partial_sum(recv_sz.begin(), recv_sz.end() - 1, recv_disp.begin() + 1);

      seq.resize(pos + S);

      MPI_Alltoallv(&send_buf[0], &tasks_sz[0], &tasks_disp[0], Type,
		    &seq[pos], &recv_sz[0], &recv_disp[0], Type, Comm);
  } // partition_balance

  template <typename Sequence>
  void partition_balance(Sequence& seq, MPI_Datatype Type, MPI_Comm Comm) {
      typedef typename Sequence::value_type value_type;
      typedef typename Sequence::iterator iterator;
      return partition_balance(seq, std::equal_to<value_type>(), detail::linear<iterator>, Type, 0, Comm);
  } // partition_balance

} // namespace mpix

#endif // MPIX2_PARTITION_BALANCE_HPP
