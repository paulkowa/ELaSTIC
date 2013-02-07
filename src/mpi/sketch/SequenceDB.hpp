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

#include <inttypes.h>

#include <mpi.h>


inline unsigned int get_nloc(int size, unsigned int n) {
    unsigned int nloc = ((static_cast<double>(n) / size) + 0.5);
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
	    MPI_Free_mem(index_);
	    MPI_Free_mem(seqs_);
	    MPI_Win_free(&index_win_);
	    MPI_Win_free(&seqs_win_);
	}
    } // ~SequenceRMA


    std::pair<bool, std::string> init(const SequenceList& SL) {
	if (is_active_ == true) return std::make_pair(false, "RMA already activated");

	int size;
	MPI_Comm_size(comm_, &size);

	i2r_ = id2rank(SL.N, size);
	unsigned int nloc = SL.seqs.size();

	// index
	if (MPI_Alloc_mem((nloc + 1) * sizeof(unsigned int), MPI_INFO_NULL, &index_) != MPI_SUCCESS) {
	    return std::make_pair(false, "could not allocate memory");
	}

	index_[0] = 0;
	for (unsigned int i = 0; i < nloc; ++i) index_[i + 1] = index_[i] + SL.seqs[i].s.size();

	MPI_Win_create(index_, (nloc + 1) * sizeof(unsigned int), sizeof(unsigned int),
		       MPI_INFO_NULL, comm_, &index_win_);

	// sequences
	if (MPI_Alloc_mem(index_[nloc], MPI_INFO_NULL, &seqs_) != MPI_SUCCESS) {
	    return std::make_pair(false, "could not allocate memory");
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

	unsigned int rank = 0;
	unsigned int offset = 0;

	boost::tie(rank, offset) = i2r_[id];

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


    template <typename Iter>
    void get(Iter first, Iter last, std::vector<std::string>& sv) {
	if (is_active_ == false) {
	    sv.clear();
	    return;
	}

	if (first == last) return;

	unsigned int rank = first->first;

	// that many sequences we need
	unsigned int m = last - first;
	last--;

	// this is how many positions they span
	unsigned int l = (last->second - first->second) + 1;

	// get index
	soffset_.resize(l + 1);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, index_win_);
	MPI_Get(&soffset_[0], l + 1, MPI_UNSIGNED, rank, first->second, l + 1, MPI_UNSIGNED, index_win_);
	MPI_Win_unlock(rank, index_win_);

	// get data
	l = soffset_.back() - soffset_.front();
	buf_.resize(l);

	MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, seqs_win_);
	MPI_Get(&buf_[0], l, MPI_CHAR, rank, soffset_[0], l, MPI_CHAR, seqs_win_);
	MPI_Win_unlock(rank, seqs_win_);

	// extract sequences
	unsigned int off = soffset_[0];
	for (unsigned int i = 0; i < soffset_.size(); ++i) soffset_[i] -= off;

	sv.resize(m);

	for (unsigned int i = 0; i < m; ++i) {
	    unsigned int pos = first[i].second - first[0].second;
	    sv[i] = std::string(buf_.begin() + soffset_[pos], buf_.begin() + soffset_[pos + 1]);
	}
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

inline bool sketch_compare(const sketch_id& lhs, const sketch_id& rhs) {
    return (lhs.sketch < rhs.sketch);
} // sketch_compare

inline sketch_id make_sketch_id(uint64_t sketch, unsigned int id, unsigned short int size) {
    sketch_id tmp;
    tmp.sketch = sketch;
    tmp.id = id;
    tmp.size = size;
    return tmp;
} // make_sketch_id

inline unsigned int hash_sketch_id(const sketch_id& si) {
    const char* ptr = reinterpret_cast<const char*>(&si);
    return *reinterpret_cast<const unsigned int*>(ptr + 2);
} // hash_sketch_id



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
    return tmp;
} // make_read_pair



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


inline unsigned int hash_read_pair0(const read_pair& rp) { return rp.id0; }


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
	return ((rp.id0 % 2 == 1) ? i2r_(rp.id0) : i2r_(rp.id1));
    } // operator()

private:
    id2rank i2r_;

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



class alignment_identity {
public:
    explicit alignment_identity(int m = 0, int s = 0, int g = 0, int h = 0)
	: align_(m, s, g, h) { }

    unsigned short int operator()(const std::string& s1, const std::string& s2) {
	unsigned int score;
	unsigned int length;
	unsigned int matches;
	boost::tie(score, length, matches) = align_(s1, s2);
	return (100 * matches) / score;
    } // operator

private:
    bio::global_alignment align_;

}; // class alignment_identity


class kmer_identity {
public:
    explicit kmer_identity(unsigned int k = 0, bool is_dna = true) : kf_(k, is_dna) { }

    unsigned short int operator()(const std::string& s1, const std::string& s2) {
	unsigned int S;
	unsigned int l0;
	unsigned int l1;
	boost::tie(S, l0, l1) = kf_(s1, s2);
	return (100 * S) / std::min(l0, l1);
    } // operator

private:
    bio::kmer_fraction kf_;

}; // class kmer_identity


template <typename Compare> class sequence_identity {
public:
    sequence_identity(const SequenceList& SL, unsigned int offset, Compare compare)
	: SL_(SL), offset_(offset), compare_(compare) {

    } // sequence_identity

    read_pair operator()(read_pair rp) {
	unsigned int id0 = rp.id0 - offset_;
	unsigned int id1 = rp.id1 - offset_;
	rp.count = compare_(SL_.seqs[id0].s, SL_.seqs[id1].s);
	return rp;
    } // operator()

    // case where id0 is local
    read_pair operator()(read_pair rp, const std::string& s) {
	unsigned int id0 = rp.id0 - offset_;
	rp.count = compare_(SL_.seqs[id0].s, s);
	return rp;
    } // operator()

    // case where id1 is local
    read_pair operator()(const std::string& s, read_pair rp) {
	unsigned int id1 = rp.id1 - offset_;
	rp.count = compare_(s, SL_.seqs[id1].s);
	return rp;
    } // operator()

private:
    const SequenceList& SL_;
    unsigned int offset_;

    Compare compare_;

}; // class sequence_identity

#endif // SEQUENCE_DB_HPP
