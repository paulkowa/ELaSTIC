/***
 *  $Id$
 **
 *  File: SequenceDB.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
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

#include <mpi.h>

#include "sketch_id.hpp"


#ifdef WITH_MPE
#include <mpix2/MPE_Log.hpp>
#endif // WITH_MPE


inline unsigned int get_nloc(int size, unsigned int n) {
    unsigned int nloc = static_cast<unsigned int>((static_cast<double>(n) / size) + 0.5);
    if (n <= nloc * (size - 1)) nloc = n / size;
    return nloc;
} // get_nloc

template <typename T>
inline std::pair<unsigned int, unsigned int> block(int size, int rank, unsigned int n) {
    unsigned int nloc = get_nloc(size, n);
    unsigned int offset = rank * nloc * sizeof(T);
    if (rank == size - 1) nloc = n - (rank * nloc);
    return std::make_pair(nloc, offset);
} // block

class id2rank {
public:
    explicit id2rank(unsigned int n = 1, int size = 1) : last_(size - 1) {
	nloc_ = get_nloc(size, n);
    } // id2rank

    // return host processor for given sequence
    unsigned int operator()(unsigned int id) const { return std::min<unsigned int>(id / nloc_, last_); }

    // return host and local id for given sequence
    std::pair<unsigned int, unsigned int> operator[](unsigned int id) const {
	unsigned int rank = this->operator()(id);
	return std::make_pair(rank, id - rank * nloc_);
    } // operator()

private:
    int last_;
    unsigned int nloc_;

}; // class id2rank



struct Sequence {
    unsigned int id;
    std::string s;
}; // struct Sequence

inline bool operator<(const Sequence& s0, const Sequence& s1) {
    return (s0.id < s1.id);
} // operator<


struct SequenceList {
    unsigned int N; // global size
    std::vector<Sequence> seqs;
}; // struct SequenceList



class SequenceRMA {
public:
    SequenceRMA(MPI_Comm comm) : comm_(comm), is_active_(false), index_(0), seqs_(0) { }

    ~SequenceRMA() {
	if (is_active_ == true) {
	    MPI_Win_free(&index_win_);
	    MPI_Win_free(&seqs_win_);
	    MPI_Free_mem(index_);
	    MPI_Free_mem(seqs_);
	}
    } // ~SequenceRMA

    std::pair<bool, std::string> init(const SequenceList& SL) {
#ifdef WITH_MPE
	mpix::MPE_Log mpe_log("SequenceRMA init", "white");
	mpe_log_.init("SequenceRMA get", "white");
#endif // WITH_MPE

	if (is_active_ == true) return std::make_pair(false, "SequenceRMA already activated");

	int size;
	MPI_Comm_size(comm_, &size);

	i2r_ = id2rank(SL.N, size);
	unsigned int nloc = SL.seqs.size();

	// index
	if (MPI_Alloc_mem((nloc + 1) * sizeof(unsigned int), MPI_INFO_NULL, &index_) != MPI_SUCCESS) {
	    return std::make_pair(false, "SequenceRMA could not allocate memory");
	}

	index_[0] = 0;
	for (unsigned int i = 0; i < nloc; ++i) index_[i + 1] = index_[i] + SL.seqs[i].s.size();

	MPI_Win_create(index_, (nloc + 1) * sizeof(unsigned int), sizeof(unsigned int),
		       MPI_INFO_NULL, comm_, &index_win_);

	// sequences
	if (MPI_Alloc_mem(index_[nloc], MPI_INFO_NULL, &seqs_) != MPI_SUCCESS) {
	    return std::make_pair(false, "SequenceRMA could not allocate memory");
	}

	char* pos = seqs_;

	for (unsigned int i = 0; i < nloc; ++i) {
	    unsigned int l = SL.seqs[i].s.size();
	    std::memcpy(pos, SL.seqs[i].s.c_str(), l);
	    pos += l;
	}

	MPI_Win_create(seqs_, index_[nloc], 1, MPI_INFO_NULL, comm_, &seqs_win_);

	is_active_ = true;
	return std::make_pair(true, "");
    } // init


    std::string get(unsigned int id) {
	if (is_active_ == false) return "";

#ifdef WITH_MPE
	mpe_log_.start();
#endif // WITH_MPE

	unsigned int rank = 0;
	unsigned int offset = 0;

	boost::tie(rank, offset) = i2r_[id];

	unsigned int len[2];

	MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, index_win_);
	MPI_Get(len, 2, MPI_UNSIGNED, rank, offset, 2, MPI_UNSIGNED, index_win_);
	MPI_Win_unlock(rank, index_win_);

	unsigned int l = len[1] - len[0];
	buf_.resize(l);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, seqs_win_);
	MPI_Get(&buf_[0], l, MPI_CHAR, rank, len[0], l, MPI_CHAR, seqs_win_);
	MPI_Win_unlock(rank, seqs_win_);

#ifdef WITH_MPE
	mpe_log_.stop();
#endif // WITH_MPE

	return std::string(buf_.begin(), buf_.end());
    } // get


    std::string get(unsigned int rank, unsigned int offset) {
	if (is_active_ == false) return "";

#ifdef WITH_MPE
	mpe_log_.start();
#endif // WITH_MPE

	unsigned int len[2];

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, index_win_);
	MPI_Get(len, 2, MPI_UNSIGNED, rank, offset, 2, MPI_UNSIGNED, index_win_);
	MPI_Win_unlock(rank, index_win_);

	unsigned int l = len[1] - len[0];
	buf_.resize(l);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, seqs_win_);
	MPI_Get(&buf_[0], l, MPI_CHAR, rank, len[0], l, MPI_CHAR, seqs_win_);
	MPI_Win_unlock(rank, seqs_win_);

#ifdef WITH_MPE
	mpe_log_.stop();
#endif // WITH_MPE

	return std::string(buf_.begin(), buf_.end());
    } // get


    template <typename Iter>
    void get(Iter first, Iter last, std::vector<std::string>& sv) {
	sv.clear();

	if (is_active_ == false) return;
	if (first == last) return;

#ifdef WITH_MPE
	mpe_log_.start();
#endif // WITH_MPE

	unsigned int rank = first->first;
	last--;

	// this is how many sequences we need
	unsigned int l = last->second - first->second + 1;

	// get index
	soffset_.resize(l + 1);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, index_win_);
	MPI_Get(&soffset_[0], l + 1, MPI_UNSIGNED, rank, first->second, l + 1, MPI_UNSIGNED, index_win_);
	MPI_Win_unlock(rank, index_win_);

	// get data
	unsigned int m = soffset_.back() - soffset_.front();
	buf_.resize(m);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, seqs_win_);
	MPI_Get(&buf_[0], m, MPI_CHAR, rank, soffset_[0], m, MPI_CHAR, seqs_win_);
	MPI_Win_unlock(rank, seqs_win_);

	// extract sequences
	unsigned int off = soffset_[0];
	for (unsigned int i = 0; i < l + 1; ++i) soffset_[i] -= off;

	sv.resize(l);

	for (unsigned int i = 0; i < l; ++i) {
	    sv[i] = std::string(buf_.begin() + soffset_[i], buf_.begin() + soffset_[i + 1]);
	}

#ifdef WITH_MPE
	mpe_log_.stop();
#endif // WITH_MPE
    } // get


private:
    SequenceRMA(const SequenceRMA&);
    void operator=(const SequenceRMA&);

    id2rank i2r_;

    MPI_Comm comm_;
    bool is_active_;

    unsigned int* index_;
    char* seqs_;

    MPI_Win index_win_;
    MPI_Win seqs_win_;

    std::vector<unsigned int> soffset_;
    std::vector<char> buf_;

#ifdef WITH_MPE
    mpix::MPE_Log mpe_log_;
#endif // WITH_MPE

}; // SequenceRMA


typedef std::vector<uint64_t> shingle_list_type;


struct read_pair {
    unsigned int id0;
    unsigned int id1;
    unsigned short int size;
    unsigned short int count;
    int score;
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

inline bool source_compare(const read_pair& lhs, const read_pair& rhs) {
    return (lhs.id0 < rhs.id0);
} // operator<

inline read_pair make_read_pair(const sketch_id& r0, const sketch_id& r1) {
    read_pair tmp;
    tmp.id0 = r0.id;
    tmp.id1 = r1.id;
    if (tmp.id1 < tmp.id0) std::swap(tmp.id0, tmp.id1);
    tmp.size = std::max<unsigned short int>(std::min(r0.size, r1.size), 1);
    tmp.count = 1;
    tmp.score = 0;
    return tmp;
} // make_read_pair

std::ostream& write_read_pair(std::ostream& os, const read_pair& rp) {
    os << rp.id0 << " " << rp.id1 << "\t" << rp.score << " " << rp.size << " " << rp.count;
    return os;
} // write_read_pair


typedef std::pair<unsigned int, unsigned int> rank_read_t;


class read2rank {
public:
    explicit read2rank(int rank, const id2rank& i2r) : rank_(rank), i2r_(i2r) { }

    rank_read_t operator()(const read_pair& rp) const {
	unsigned int rank0, id0;
	unsigned int rank1, id1;
	boost::tie(rank0, id0) = i2r_[rp.id0];
	boost::tie(rank1, id1) = i2r_[rp.id1];
	return (rank_ == rank0) ? std::make_pair(rank1, id1) : std::make_pair(rank0, id0);
    } // operator

private:
    int rank_;
    const id2rank& i2r_;

}; // read2rank


class compare_rank {
public:
    explicit compare_rank(const read2rank& r2r) : r2r_(r2r) { }

    bool operator()(const read_pair& rd0, const read_pair& rd1) const {
	unsigned int rank0, id0;
	unsigned int rank1, id1;
	boost::tie(rank0, id0) = r2r_(rd0);
	boost::tie(rank1, id1) = r2r_(rd1);
	return ((rank0 < rank1) || (!(rank1 < rank0) && (id0 < id1)));
    } // operator()

private:
    read2rank r2r_;

}; // class compare_rank


inline unsigned int hash_read_pair0(const read_pair& rp) { return rp.id0 ^ rp.id1; }


class hash_read_pair1 {
public:
    hash_read_pair1(unsigned int n, int size) : i2r_(n, size) { }

    unsigned int operator()(const read_pair& rp) const { return i2r_(rp.id0); }

private:
    id2rank i2r_;

}; // hash_read_pair1


class hash_read_pair2 {
public:
    hash_read_pair2(unsigned int n, int size) : i2r_(n, size) { }

    unsigned int operator()(const read_pair& rp) const {
	return (rp.id0 + rp.id1) % 2 ? i2r_(rp.id0) : i2r_(rp.id1);
    } // operator()

private:
    id2rank i2r_;

}; // hash_read_pair2


class local {
public:
    local(unsigned int lo, unsigned int hi) : lo_(lo), hi_(hi) { }

    bool operator()(const read_pair& rp) const {
	return !((rp.id0 < lo_) || (rp.id0 > hi_) || (rp.id1 < lo_) || (rp.id1 > hi_));
    } // operator()

private:
    unsigned int lo_;
    unsigned int hi_;

}; // class local


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
    not_collocated(unsigned int n, int size) :  i2r_(n, size) { }

    bool operator()(const read_pair& rp) const {
	int p0 = i2r_(rp.id0);
	int p1 = i2r_(rp.id1);
	return (p0 != p1);
    } // operator

private:
    id2rank i2r_;

}; // class not_collocated


inline read_pair kmer_fraction(read_pair rp) {
    rp.count = (100 * rp.count) / rp.size;
    return rp;
} // kmer_fraction


class appx_kmer_fraction {
public:
    explicit appx_kmer_fraction(unsigned short int kmer, short int mod)
	: kmer_(kmer), mod_(mod) { }

    // this works when size_ is actual sequence length
    read_pair operator()(read_pair rp) const {
	rp.count = (100 * mod_ * rp.count) / (rp.size - kmer_ + 1);
	return rp;
    } // operator()

private:
    unsigned short int kmer_;
    short int mod_;

}; // class appx_kmer_fraction


inline unsigned short int get_kmer_fraction(const read_pair& rp) {
    return (100 * rp.score) / std::min(rp.size, rp.count);
} // get_kmer_fraction

inline unsigned short int get_alignment_identity(const read_pair& rp) {
    // this is bit messy :-)
    // for cd-hit identity rp.size is corrected by compare
    // (called inside sequence_compare)
    return (100 * rp.count) / rp.size;
} // get_alignment_identity


class not_similar {
public:
    explicit not_similar(unsigned short int jmin) : jmin_(jmin) { }

    bool operator()(const read_pair& rp) const {
        return rp.count < jmin_;
    } // operator()

private:
    unsigned short int jmin_;

}; // class not_similar


class not_similar2 {
public:
    not_similar2(unsigned int method, unsigned short int jmin)
	: method_(method), jmin_(jmin) { }

    bool operator()(const read_pair& rp) const {
	unsigned int score = 0;
	if (method_ == 0) {
	    score = get_kmer_fraction(rp);
	} else {
	    score = get_alignment_identity(rp);
	}

	return score < jmin_;
    } // operator()

private:
    unsigned int method_;
    unsigned short int jmin_;

}; // class not_similar2


class read_pair_count {
public:
    explicit read_pair_count(unsigned int method) : method_(method) { }

    read_pair operator()(read_pair rp) const {
	unsigned int score = 0;
	if (method_ == 0) {
	    score = get_kmer_fraction(rp);
	} else {
	    score = get_alignment_identity(rp);
	}
	rp.count = score;
	return rp;
    } // operator()

private:
    unsigned int method_;

}; // class read_pair_count


template <typename Compare> class sequence_compare {
public:
    sequence_compare(const SequenceList& SL, unsigned int offset, Compare compare)
	: SL_(SL), offset_(offset), compare_(compare) { }

    read_pair operator()(read_pair rp) {
	unsigned int id0 = rp.id0 - offset_;
	unsigned int id1 = rp.id1 - offset_;
	boost::tuple<int, int, int> cmp = compare_(SL_.seqs[id0].s, SL_.seqs[id1].s);
	prv_assign__(cmp, rp);
	return rp;
    } // operator()

    // case where id0 is local
    read_pair operator()(read_pair rp, const std::string& s) {
	unsigned int id0 = rp.id0 - offset_;
	boost::tuple<int, int, int> cmp = compare_(SL_.seqs[id0].s, s);
	prv_assign__(cmp, rp);
	return rp;
    } // operator()

    // case where id1 is local
    read_pair operator()(const std::string& s, read_pair rp) {
	unsigned int id1 = rp.id1 - offset_;
	boost::tuple<int, int, int> cmp = compare_(s, SL_.seqs[id1].s);
	prv_assign__(cmp, rp);
	return rp;
    } // operator()

    // case where none is local
    read_pair operator()(read_pair rp, const std::string& s0, const std::string& s1) {
	boost::tuple<int, int, int> cmp = compare_(s0, s1);
	prv_assign__(cmp, rp);
	return rp;
    } // operator()

private:
    void prv_assign__(const boost::tuple<int, int, int>& cmp, read_pair& rp) {
	rp.score = boost::get<0>(cmp);
	rp.size =  boost::get<1>(cmp);
	rp.count = boost::get<2>(cmp);
    } // prv_assign__

    const SequenceList& SL_;
    unsigned int offset_;

    Compare compare_;

}; // class sequence_identity

#endif // SEQUENCE_DB_HPP
