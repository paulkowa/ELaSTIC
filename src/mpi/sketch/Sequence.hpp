/***
 *  $Id$
 **
 *  File: Sequence.hpp
 *  Created: Dec 16, 2014
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2014 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <vector>


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


typedef std::vector<uint64_t> shingle_list_type;

#endif // SEQUENCE_HPP
