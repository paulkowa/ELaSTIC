/***
 *  $Id$
 **
 *  File: Sequence.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <iostream>
#include <string>
#include <vector>

#include <bio/sequence_compare.hpp>
#include <jaz/hash.hpp>

#include <inttypes.h>


struct Sequence {
    unsigned int id;
    std::string s;
}; // struct Sequence

inline bool operator<(const Sequence& s0, const Sequence& s1) {
    return (s0.id < s1.id);
} // operator<


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


class approx_jaccard {
public:
    explicit approx_jaccard(unsigned short int kmer, short int mod)
	: kmer_(kmer), mod_(mod) { }

    read_pair operator()(read_pair rp) const {
	rp.count = (100 * mod_ * rp.count) / (rp.size - kmer_ + 1);
	return rp;
    } // operator()

private:
    unsigned short int kmer_;
    short int mod_;

}; // class approx_jaccard


class identity {
public:
    identity(const std::vector<Sequence>& seqs,
	     unsigned short int method,
	     unsigned short int kmer)
	: seqs_(seqs), method_(method), A_(5, -4, -10, -1), J_(kmer) { }

    read_pair operator()(read_pair rp) {
	s0_.id = rp.id0;
	s1_.id = rp.id1;

	std::vector<Sequence>::const_iterator id0, id1;
	id0 = std::equal_range(seqs_.begin(), seqs_.end(), s0_).first;
	id1 = std::equal_range(seqs_.begin(), seqs_.end(), s1_).first;

	if (method_ == 0) rp.count = A_(id0->s, id1->s);
	else rp.count = J_(id0->s, id1->s);

	return rp;
    } // operator()

private:
    const std::vector<Sequence>& seqs_;
    unsigned short int method_;

    bio::global_alignment A_;
    bio::dna_jaccard_index J_;

    Sequence s0_;
    Sequence s1_;

}; // class identity


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

}; // class is_local


class not_collocated {
public:
    not_collocated(unsigned int n, int size) : nloc_(n / size), last_(size - 1) { }

    bool operator()(const read_pair& rp) const {
	unsigned int p0 = std::min((rp.id0 / nloc_), last_);
	unsigned int p1 = std::min((rp.id1 / nloc_), last_);
	return (p0 != p1);
    } // operator

private:
    unsigned int nloc_;
    unsigned int last_;

}; // class not_collocated


inline read_pair make_read_pair(const sketch_id& r0, const sketch_id& r1) {
    read_pair tmp;
    tmp.id0 = r0.id;
    tmp.id1 = r1.id;
    if (tmp.id1 < tmp.id0) std::swap(tmp.id0, tmp.id1);
    tmp.size = std::min(r0.size, r1.size);
    tmp.count = 1;
    return tmp;
} // make_read_pair


inline unsigned int read_pair_hash(const read_pair& rp) { return rp.id0; }

class read_pair_local_hash {
public:
    read_pair_local_hash(unsigned int n, int size) : nloc_(n / size), last_(size - 1) { }

    unsigned int operator()(const read_pair& rp) const {
	unsigned int p = std::min((rp.id0 / nloc_), last_);
	return p;
    } // operator

private:
    unsigned int nloc_;
    unsigned int last_;

}; // class not_collocated

#endif // SEQUENCE_HPP
