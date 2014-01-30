/***
 *  $Id$
 **
 *  File: sketch_id.hpp
 *  Created: Oct 15, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef SKETCH_ID_HPP
#define SKETCH_ID_HPP

#include <algorithm>
#include <iostream>
#include <vector>

#include <jaz/algorithm.hpp>

#include <inttypes.h>


struct sketch_id {
    uint64_t sketch;
    unsigned int id;
    unsigned short int size;
    unsigned short int sep;
    unsigned int part;
}; // struct sketch_id

inline std::ostream& operator<<(std::ostream& os, const sketch_id& si) {
    os << si.sketch << " " << si.id << " " << si.size << " " << si.sep << " " << si.part;
    return os;
} // operator<<

inline bool operator==(const sketch_id& lhs, const sketch_id& rhs) {
    return (lhs.sketch == rhs.sketch) && (lhs.part == rhs.part);
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
    tmp.sep = 0;
    tmp.part = 0;
    return tmp;
} // make_sketch_id

inline unsigned int hash_sketch_id(const sketch_id& si) {
    const char* ptr = reinterpret_cast<const char*>(&si);
    return *reinterpret_cast<const unsigned int*>(ptr + 2);
} // hash_sketch_id

inline unsigned int hash_sketch_id2(const sketch_id& si) {
    uint64_t x = si.sketch;
    x >>= 3;
    return (x ^ (x >> 10) ^ (x >> 20));
} // hash_sketch_id2

struct id_sketch {
    unsigned int id;
    uint64_t sketch;
}; // struct id_sketch

inline bool operator<(const id_sketch& lhs, const id_sketch& rhs) {
    return ((lhs.id < rhs.id) || (!(rhs.id < lhs.id) && (lhs.sketch < rhs.sketch)));
} // operator<

inline bool id_compare2(const id_sketch& lhs, const id_sketch& rhs) {
    return (lhs.id < rhs.id);
} // id_compare

inline bool sketch_compare2(const id_sketch& lhs, const id_sketch& rhs) {
    return (lhs.sketch < rhs.sketch);
} // sketch_compare

inline id_sketch make_id_sketch(unsigned int id, uint64_t sketch) {
    id_sketch tmp;
    tmp.id = id;
    tmp.sketch = sketch;
    return tmp;
} // make_id_sketch


template <typename Iter>
inline int sketch_part_cost(Iter first, Iter last) {
    int n = (last - first);
    if (first->sep == 0) return (n * (n - 1)) >> 1;
    if (first->sep == n) return n * n;
    return first->sep * (n - first->sep);
} // sketch_part_cost


inline void decompose_sketch_triangle(std::vector<sketch_id>& sketch_list, unsigned int& part, int l) {
    int pos = sketch_list.size() - l;
    int p = l / 2;

    // top left stays unchanged

    // bottom right becomes new partition
    for (int i = pos + p; i != pos + l; ++i) sketch_list[i].part = part;
    part++;

    // top right rectangle
    for (int i = pos; i != pos + l; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = p;
    }
    part++;
} // decompose_sketch_triangle

inline void decompose_sketch_rectangle(std::vector<sketch_id>& sketch_list, unsigned int& part, int l) {
    int pos = sketch_list.size() - l;

    int p = sketch_list[pos].sep;
    int px = p / 2;
    int py = (l - p) / 2;

    // top left
    for (int i = pos; i != pos + px; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().sep = px;
    }
    for (int i = pos + p; i != pos + p + py; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().sep = px;
    }
    // part++;

    // top right
    for (int i = pos; i != pos + px; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = px;
    }
    for (int i = pos + p + py; i != pos + l; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = px;
    }
    part++;

    // bottom left
    for (int i = pos + px; i != pos + p; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = p - px;
    }
    for (int i = pos + p; i != pos + p + py; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = p - px;
    }
    part++;

    // bottom right
    for (int i = pos + px; i != pos + p; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = p - px;
    }
    for (int i = pos + p + py; i != pos + l; ++i) {
	sketch_list.push_back(sketch_list[i]);
	sketch_list.back().part = part;
	sketch_list.back().sep = p - px;
    }
    part++;

    std::rotate(sketch_list.begin() + pos, sketch_list.begin() + pos + l, sketch_list.end());
    sketch_list.resize(sketch_list.size() - l);
} // decompose_sketch_rectangle

inline unsigned int decompose_sketch_list(std::vector<sketch_id>& sketch_list, int lim, unsigned int part) {
    int first = 0;

    while ((sketch_list.begin() + first) != sketch_list.end()) {
	int temp = jaz::range(sketch_list.begin() + first, sketch_list.end()) - sketch_list.begin();
	int l = temp - first;

	if (lim < l) {
	    bool is_triangle = (sketch_list[first].sep == 0);

	    // get partition to the end
	    std::rotate(sketch_list.begin() + first, sketch_list.begin() + temp, sketch_list.end());

	    if (is_triangle) decompose_sketch_triangle(sketch_list, part, l);
	    else decompose_sketch_rectangle(sketch_list, part, l);

	} else first = temp;
    } // while

    return part;
} // decompose_sketch_list

#endif // SKETCH_ID_HPP
