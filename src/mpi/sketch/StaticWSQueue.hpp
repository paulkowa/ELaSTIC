/***
 *  $Id$
 **
 *  File: StaticWSQueue.hpp
 *  Created: Feb 22, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef STATIC_WS_QUEUE_HPP
#define STATIC_WS_QUEUE_HPP

#include <algorithm>
#include <mpi.h>


template <typename T> class StaticWSQueue {
public:
    struct range_type {
	unsigned int first;
	unsigned int last;
    }; // struct range_type


    static range_type make_range(unsigned int first, unsigned int last) {
	range_type range;
	range.first = first;
	range.last = last;
	return range;
    } // make_range


    explicit StaticWSQueue(MPI_Comm comm) : active_(false), comm_(comm) {
	MPI_Comm_size(comm_, &size_);
	MPI_Comm_rank(comm_, &rank_);
    } // StaticWSQueue


    ~StaticWSQueue() {
	finalize();
    } // ~StaticWSQueue


    template <typename IterData, typename IterTask>
    bool init(IterData dfirst, IterData dlast, IterTask tfirst, IterTask tlast) {
	if (active_ == true) return false;
	active_ = true;

	victims_.resize(size_);
	for (unsigned int i = 0; i < size_; ++i) victims_[i] = i;
	victims_.erase(victims_.begin() + rank_);

	ends_ = 1;

	MPI_Type_contiguous(sizeof(range_type), MPI_BYTE, &MPI_RANGE_TYPE);
	MPI_Type_commit(&MPI_RANGE_TYPE);

	unsigned int n = dlast - dfirst;

	if (MPI_Alloc_mem(n * sizeof(T), MPI_INFO_NULL, &data_) != MPI_SUCCESS) return false;
	std::copy(dfirst, dlast, data_);
	MPI_Win_create(data_, n * sizeof(T), sizeof(T), MPI_INFO_NULL, comm_, &data_win_);

	unsigned int m = tlast - tfirst;

	task_.resize(m);
	std::copy(tfirst, tlast, task_.begin());

	queue_[0] = 0;
	queue_[1] = m;

	MPI_Irecv(&req_buf_, 1, MPI_INT, MPI_ANY_SOURCE, SWSQ_STEAL_REQ, comm_, &req_);

	return true;
    } // init

    void finalize() {
	if (active_ == true) {
	    MPI_Status stat;

	    while (progress() == true) { }
	    MPI_Cancel(&req_);

	    MPI_Win_free(&data_win_);
	    MPI_Free_mem(data_);
	    MPI_Type_free(&MPI_RANGE_TYPE);
	} // if

	active_ = false;
    } // finalize


    bool progress() {
	if (size_ <= ends_) return false;

	MPI_Status stat;
	int flag = 0;

	MPI_Test(&req_, &flag, &stat);

	if (flag == true) {
	    range_type task;

	    if (queue_[1] > 0) queue_[1]--;

	    if ((queue_[1] > 0) && (queue_[0] <= queue_[1])) task = task_[queue_[1]];
	    else {
		task.first = task.last = 0;
		ends_++;
	    }

	    MPI_Send(&task, 1, MPI_RANGE_TYPE, req_buf_, SWSQ_TASK_ANS, comm_);
	    MPI_Irecv(&req_buf_, 1, MPI_INT, MPI_ANY_SOURCE, SWSQ_STEAL_REQ, comm_, &req_);
	} // if

	return (ends_ < size_);
    } // progress


    bool get(const T*& first, const T*& last) {
	first = 0;
	last = 0;

	if (active_ == false) return false;
	bool res = true;

	if (queue_[1] <= queue_[0]) res = false;
	else queue_[0]++;

	first = data_ + task_[queue_[0] - 1].first;
	last = data_ + task_[queue_[0] - 1].last;

	return res;
    } // get


    bool steal(T*& first, T*& last) {
	if (size_ == 1) return false;

	MPI_Status stat;

	MPI_Request sreq;
	MPI_Request rreq;

	int sflag = 0;
	int rflag = 0;

	int vpos = 0;
	int vic = 0;

	range_type task;

	do {
	    vpos = std::rand() % victims_.size();
	    vic = victims_[vpos];

	    task.first = task.last = 0;

	    MPI_Isend(&rank_, 1, MPI_INT, vic, SWSQ_STEAL_REQ, comm_, &sreq);
	    MPI_Irecv(&task, 1, MPI_RANGE_TYPE, vic, SWSQ_TASK_ANS, comm_, &rreq);

	    sflag = 0;

	    do {
		MPI_Test(&sreq, &sflag, &stat);
		progress();
	    } while (sflag == false);

	    rflag = 0;

	    do {
		MPI_Test(&rreq, &rflag, &stat);
		progress();
	    } while (rflag == false);

	    if (task.first == task.last) victims_.erase(victims_.begin() + vpos);
	    else {
		unsigned int k = task.last - task.first;

		first = new T[k];
		last = first + k;

		MPI_Win_lock(MPI_LOCK_SHARED, vic, MPI_MODE_NOCHECK, data_win_);
		MPI_Get(first, k * sizeof(T), MPI_BYTE, vic, task.first, k * sizeof(T), MPI_BYTE, data_win_);
		MPI_Win_unlock(vic, data_win_);
	    }
	} while ((task.first == task.last) && (victims_.empty() == false));

	if (victims_.empty() == true) return false;
	return true;
    } // steal


private:
    enum { SWSQ_STEAL_REQ = 111, SWSQ_TASK_ANS = 222 };

    bool active_;

    MPI_Comm comm_;
    int rank_;
    int size_;

    MPI_Request req_;
    int req_buf_;

    MPI_Datatype MPI_RANGE_TYPE;

    // data and task storage
    MPI_Win data_win_;
    T* data_;

    std::vector<range_type> task_;

    // queue pointers [head,tail]
    int queue_[2];

    std::vector<int> victims_;
    int ends_;

}; // class StaticWSQueue

#endif // STATIC_WS_QUEUE_HPP
