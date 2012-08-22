/***
 *  $Id$
 **
 *  File: SequenceDB.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef SEQUENCE_DB_HPP
#define SEQUENCE_DB_HPP

#include <iostream>
#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include <bio/sequence_compare.hpp>
#include <jaz/hash.hpp>

#include <inttypes.h>

#include <mpi.h>


template <typename T>
inline std::pair<unsigned int, unsigned int> block(int size, int rank, unsigned int n) {
    unsigned int nloc = ((static_cast<double>(n) / size) + 0.5);
    unsigned int offset = rank * nloc * sizeof(T);
    if (rank == size - 1) nloc = n - (rank * nloc);
    return std::make_pair(nloc, offset);
} // block

class id2cpu {
public:
    explicit id2cpu(unsigned int n = 1, int size = 1) {
	nloc_ = ((static_cast<double>(n) / size) + 0.5);
    } // id2cpu

    unsigned int operator()(unsigned int id) const { return id / nloc_; }

    std::pair<unsigned int, unsigned int> operator[](unsigned int id) const {
	unsigned int rank = this->operator()(id);
	return std::make_pair(rank, id - rank * nloc_);
    } // operator()

private:
    unsigned int nloc_;

}; // class id2cpu



struct Sequence {
    unsigned int id;
    std::string s;
}; // struct Sequence

inline bool operator<(const Sequence& s0, const Sequence& s1) {
    return (s0.id < s1.id);
} // operator<



class SequenceRMA {
public:
    SequenceRMA(MPI_Comm comm, unsigned int n)
	: comm_(comm), is_active_(false), index_(0), seqs_(0) {
	int size;
	MPI_Comm_size(comm, &size);
	i2c_ = id2cpu(n, size);
    } // SequenceRMA

    ~SequenceRMA() {
	if (is_active_ == true) {
	    MPI_Free_mem(index_);
	    MPI_Free_mem(seqs_);
	    MPI_Win_free(&index_win_);
	    MPI_Win_free(&seqs_win_);
	}
    } // ~SequenceRMA


    std::pair<bool, std::string> init(const std::vector<Sequence>& seqs) {
	if (is_active_ == true) return std::make_pair(false, "RMA already activated");

	unsigned int nloc = seqs.size();

	// index
	if (MPI_Alloc_mem((nloc + 1) * sizeof(unsigned int), MPI_INFO_NULL, &index_) != MPI_SUCCESS) {
	    return std::make_pair(false, "could not allocate memory");
	}

	index_[0] = 0;
	for (unsigned int i = 0; i < nloc; ++i) index_[i + 1] = index_[i] + seqs[i].s.size();

	MPI_Win_create(index_, (nloc + 1) * sizeof(unsigned int), sizeof(unsigned int),
		       MPI_INFO_NULL, comm_, &index_win_);

	// sequences
	if (MPI_Alloc_mem(index_[nloc], MPI_INFO_NULL, &seqs_) != MPI_SUCCESS) {
	    return std::make_pair(false, "could not allocate memory");
	}

	char* pos = seqs_;

	for (unsigned int i = 0; i < nloc; ++i) {
	    unsigned int l = seqs[i].s.size();
	    std::memcpy(pos, seqs[i].s.c_str(), l);
	    pos += l;
	}

	MPI_Win_create(seqs_, index_[nloc], 1, MPI_INFO_NULL, comm_, &seqs_win_);

	is_active_ = true;
	return std::make_pair(true, "");
    } // init


    std::string get(unsigned int id) {
	if (is_active_ == false) return "";

	unsigned int rank = 0;
	unsigned int offset = 0;
	boost::tie(rank, offset) = i2c_[id];

	unsigned int len[2];

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, index_win_);
	MPI_Get(len, 2, MPI_UNSIGNED, rank, offset, 2, MPI_UNSIGNED, index_win_);
	MPI_Win_unlock(rank, index_win_);

	unsigned int l = len[1] - len[0];
	buf_.resize(l);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, seqs_win_);
	MPI_Get(&buf_[0], l, MPI_CHAR, rank, len[0], l, MPI_CHAR, seqs_win_);
	MPI_Win_unlock(rank, seqs_win_);

	return std::string(buf_.begin(), buf_.end());
    } // get


private:
    SequenceRMA(const SequenceRMA&);
    void operator=(const SequenceRMA&);

    id2cpu i2c_;

    MPI_Comm comm_;
    bool is_active_;

    unsigned int* index_;
    char* seqs_;

    MPI_Win index_win_;
    MPI_Win seqs_win_;

    std::vector<char> buf_;

}; // SequenceRMA



typedef std::vector<uint64_t> shingle_list_type;



struct sketch_id {
    uint64_t sketch;
    unsigned int id;
    unsigned short int size;
}; // struct sketch_read

inline std::ostream& operator<<(std::ostream& os, const sketch_id& si) {
    os << si.sketch << " " << si.id << " " << si.size;
    return os;
} // operator<<

inline bool operator==(const sketch_id& lhs, const sketch_id& rhs) {
    return lhs.sketch == rhs.sketch;
} // operator==

inline bool operator!=(const sketch_id& lhs, const sketch_id& rhs) {
    return !(lhs == rhs);
} // operator!=

inline bool operator<(const sketch_id& lhs, const sketch_id& rhs) {
    return ((lhs.sketch < rhs.sketch) || (!(rhs.sketch < lhs.sketch) && (lhs.id < rhs.id)));
} // operator<

inline sketch_id make_sketch_id(uint64_t sketch, unsigned int id, unsigned short int size) {
    sketch_id tmp;
    tmp.sketch = sketch;
    tmp.id = id;
    tmp.size = size;
    return tmp;
} // make_sketch_id

inline unsigned int sketch_id_hash(const sketch_id& si) {
    const char* ptr = reinterpret_cast<const char*>(&si);
    return *reinterpret_cast<const unsigned int*>(ptr + 2);
} // sketch_id_hash



struct read_pair {
    unsigned int id0;
    unsigned int id1;
    unsigned short int size;
    unsigned short int count;
}; // struct read_pair

inline std::ostream& operator<<(std::ostream& os, const read_pair& rp) {
    os << rp.id0 << " " << rp.id1 << "\t" << rp.count;
    return os;
} // operator<<

inline bool operator==(const read_pair& lhs, const read_pair& rhs) {
    return (lhs.id0 == rhs.id0) && (lhs.id1 == rhs.id1);
} // operator==

inline bool operator!=(const read_pair& lhs, const read_pair& rhs) {
    return !(lhs == rhs);
} // operator!=

inline bool operator<(const read_pair& lhs, const read_pair& rhs) {
    return ((lhs.id0 < rhs.id0) || (!(rhs.id0 < lhs.id0) && (lhs.id1 < rhs.id1)));
} // operator<

inline read_pair operator+(const read_pair& lhs, const read_pair& rhs) {
    read_pair tmp(lhs);
    tmp.count += rhs.count;
    return tmp;
} // operator+

inline read_pair make_read_pair(const sketch_id& r0, const sketch_id& r1) {
    read_pair tmp;
    tmp.id0 = r0.id;
    tmp.id1 = r1.id;
    if (tmp.id1 < tmp.id0) std::swap(tmp.id0, tmp.id1);
    tmp.size = std::min(r0.size, r1.size);
    tmp.count = 1;
    return tmp;
} // make_read_pair

inline unsigned int hash_read_pair0(const read_pair& rp) { return rp.id0; }

class hash_read_pair1 {
public:
    hash_read_pair1(unsigned int n, int size) : i2c_(n, size) { }

    unsigned int operator()(const read_pair& rp) const { return i2c_(rp.id0); }

private:
    id2cpu i2c_;

}; // hash_read_pair1

class hash_read_pair2 {
public:
    hash_read_pair2(unsigned int n, int size) : i2c_(n, size) { }

    unsigned int operator()(const read_pair& rp) const {
	return ((rp.id0 % 2 == 1) ? i2c_(rp.id0) : i2c_(rp.id1));
    } // operator()

private:
    id2cpu i2c_;

}; // hash_read_pair2


class not_similar {
public:
    explicit not_similar(unsigned short int jmin) : jmin_(jmin) { }

    bool operator()(const read_pair& rp) const {
	return rp.count < jmin_;
    } // operator()

private:
    unsigned short int jmin_;

}; // class not_similar


class not_local {
public:
    not_local(unsigned int lo, unsigned int hi) : lo_(lo), hi_(hi) { }

    bool operator()(const read_pair& rp) const {
	return ((rp.id0 < lo_) || (rp.id0 > hi_) || (rp.id1 < lo_) || (rp.id1 > hi_));
    } // operator()

private:
    unsigned int lo_;
    unsigned int hi_;

}; // class not_local


class partially_local {
public:
    partially_local(unsigned int lo, unsigned int hi) : lo_(lo), hi_(hi) { }

    bool operator()(const read_pair& rp) const {
	return (((rp.id0 >= lo_) && (rp.id0 <= hi_)) || ((rp.id1 >= lo_) && (rp.id1 <= hi_)));
    } // operator()

private:
    unsigned int lo_;
    unsigned int hi_;

}; // class partially_local


class not_collocated {
public:
    not_collocated(unsigned int n, int size) :  i2c_(n, size) { }

    bool operator()(const read_pair& rp) const {
	int p0 = i2c_(rp.id0);
	int p1 = i2c_(rp.id1);
	return (p0 != p1);
    } // operator

private:
    id2cpu i2c_;

}; // class not_collocated



class approximate_jaccard {
public:
    explicit approximate_jaccard(unsigned short int kmer, short int mod)
	: kmer_(kmer), mod_(mod) { }

    read_pair operator()(read_pair rp) const {
	rp.count = (100 * mod_ * rp.count) / (rp.size - kmer_ + 1);
	return rp;
    } // operator()

private:
    unsigned short int kmer_;
    short int mod_;

}; // class approximate_jaccard


class sequence_identity {
public:
    sequence_identity(const std::vector<Sequence>& seqs,
		      unsigned short int method, unsigned short int kmer,
		      int size, int rank, unsigned int n)
	: seqs_(seqs), method_(method), A_(5, -4, -10, -1), J_(kmer) {
	offset_ = rank * static_cast<unsigned int>((static_cast<double>(n) / size) + 0.5);
    } // sequence_identity

    read_pair operator()(read_pair rp) {
	unsigned int id0 = rp.id0 - offset_;
	unsigned int id1 = rp.id1 - offset_;

	if (method_ == 0) rp.count = A_(seqs_[id0].s, seqs_[id1].s);
	else rp.count = J_(seqs_[id0].s, seqs_[id1].s);

	return rp;
    } // operator()

    // case where id0 is local
    read_pair operator()(read_pair rp, const std::string& s) {
	unsigned int id0 = rp.id0 - offset_;

	if (method_ == 0) rp.count = A_(seqs_[id0].s, s);
	else rp.count = J_(seqs_[id0].s, s);

	return rp;
    } // operator()

    // case where id1 is local
    read_pair operator()(const std::string& s, read_pair rp) {
	unsigned int id1 = rp.id1 - offset_;

	if (method_ == 0) rp.count = A_(s, seqs_[id1].s);
	else rp.count = J_(s, seqs_[id1].s);

	return rp;
    } // operator()

private:
    const std::vector<Sequence>& seqs_;
    unsigned short int method_;

    unsigned int offset_;

    bio::global_alignment A_;
    bio::dna_jaccard_index J_;

}; // class sequence_identity

#endif // SEQUENCE_DB_HPP
