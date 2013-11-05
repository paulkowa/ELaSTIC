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


    ~StaticWSQueue() { finalize(); }


    template <typename IterData, typename IterTask>
    bool init(IterData dfirst, IterData dlast, IterTask tfirst, IterTask tlast) {
	if (active_ == true) return false;
	active_ = true;

	victims_.resize(size_);
	for (int i = 0; i < size_; ++i) victims_[i] = i;

	victims_.erase(victims_.begin() + rank_);
	ends_ = 1;

	MPI_Type_contiguous(sizeof(range_type), MPI_BYTE, &MPI_RANGE_TYPE);
	MPI_Type_commit(&MPI_RANGE_TYPE);

	int n = dlast - dfirst;

	data_.resize(n);
	std::copy(dfirst, dlast, data_.begin());

	int m = tlast - tfirst;

	task_.resize(m);
	std::copy(tfirst, tlast, task_.begin());

	queue_[0] = -1;
	queue_[1] = m;

	MPI_Irecv(req_buf_, 2, MPI_INT, MPI_ANY_SOURCE, SWSQ_STEAL_REQ, comm_, &req_);

	return true;
    } // init


    void finalize() {
	if (active_ == true) {
	    terminate();
	    while (progress() == true) { }
	    MPI_Cancel(&req_);
	} // if

	active_ = false;
    } // finalize


    bool get(const T*& first, const T*& last) {
	first = 0;
	last = 0;

	if (active_ == false) return false;

	if ((queue_[0] + 1) == queue_[1]) return false;
	else queue_[0]++;

	first = &data_[task_[queue_[0]].first];
	last = &data_[task_[queue_[0]].last];

	return true;
    } // get


    bool steal(int& vrank, T*& first, T*& last) {
	if (size_ == 1) return false;
	vrank = -1;

	MPI_Status stat;

	MPI_Request sreq;
	MPI_Request rreq;

	int sflag = 0;
	int rflag = 0;

	int vpos = 0;
	int vic = 0;

	range_type task;

	int quest[2];
	quest[0] = rank_;
	quest[1] = rank_;

	do {
	    vpos = std::rand() % victims_.size();
	    vic = victims_[vpos];

	    task.first = task.last = 0;

	    MPI_Isend(quest, 2, MPI_INT, vic, SWSQ_STEAL_REQ, comm_, &sreq);
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
		int k = task.last - task.first;

		first = new T[k];
		last = first + k;

		MPI_Recv(first, k * sizeof(T), MPI_BYTE, vic, SWSQ_TASK_DATA, comm_, &stat);
		vrank = vic;
	    }
	} while ((task.first == task.last) && (victims_.empty() == false));

	if (victims_.empty() == true) return false;
	return true;
    } // steal


    bool progress() {
	if (size_ <= ends_) return false;

	MPI_Status stat;
	int flag = 0;

	unsigned int t = 0;

	do {
	    MPI_Test(&req_, &flag, &stat);

	    if (flag == true) {
		range_type task;
		bool ht = true;

		if (req_buf_[1] != -1) {
		    if ((queue_[0] + 1) < queue_[1]) {
			queue_[1]--;
			task = task_[queue_[1]];
		    } else ht = false;
		} else ht = false;

		if (ht == false) {
		    task.first = task.last = 0;
		    ends_++;
		}

		MPI_Send(&task, 1, MPI_RANGE_TYPE, req_buf_[0], SWSQ_TASK_ANS, comm_);
		if (ht == true) {
		    unsigned int k = task.last - task.first;
		    MPI_Send(&data_[task.first], k * sizeof(T), MPI_BYTE, req_buf_[0], SWSQ_TASK_DATA, comm_);
		}
		MPI_Irecv(req_buf_, 2, MPI_INT, MPI_ANY_SOURCE, SWSQ_STEAL_REQ, comm_, &req_);

		t++;
	    } // if
	} while ((flag == true) && (t < 1));

	return (ends_ < size_);
    } // progress


    void terminate() {
	MPI_Status stat;

	MPI_Request sreq;
	MPI_Request rreq;

	int sflag = 0;
	int rflag = 0;

	range_type task;
	int quest[2];

	quest[0] = rank_;
	quest[1] = -1;

	for (unsigned int i = 0; i < victims_.size(); ++i) {
	    task.first = task.last = 0;

	    MPI_Isend(quest, 2, MPI_INT, victims_[i], SWSQ_STEAL_REQ, comm_, &sreq);
	    MPI_Irecv(&task, 1, MPI_RANGE_TYPE, victims_[i], SWSQ_TASK_ANS, comm_, &rreq);

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
	}
    } // terminate


private:
    enum { SWSQ_STEAL_REQ = 111, SWSQ_TASK_ANS = 222, SWSQ_TASK_DATA = 333 };

    bool active_;

    MPI_Comm comm_;
    int rank_;
    int size_;

    MPI_Request req_;
    int req_buf_[2];

    MPI_Datatype MPI_RANGE_TYPE;

    // data and task storage
    std::vector<T> data_;

    std::vector<range_type> task_;

    // queue pointers [head,tail]
    int queue_[2];

    std::vector<int> victims_;
    int ends_;

}; // class StaticWSQueue

#endif // STATIC_WS_QUEUE_HPP
