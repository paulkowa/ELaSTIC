/***
 *  $Id$
 **
 *  File: sffun.h
 *  Created: Aug 28, 2008
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2007-2011 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of jaz.
 */

#ifndef SFFUN_H
#define SFFUN_H

#include <limits.h>
#include <math.h>


// See http://www.research.att.com/~njas/sequences/a073189.txt
// for triangular and rectangular SF functions.

// lower-triangle matrix representation (0)
//
//    0  1  2  3  4  j
//   -- -- -- -- -- --
// 0|
// 1| 0
// 2| 1  2
// 3| 3  4  5
// 4| 6  7  8  9
// i|         n
//
// n = ((i * i) - i) / 2 + j
// i = floor(1/2 + sqrt(2n + 2))
// j = n - (i * i - i) / 2

inline unsigned int ltsf0(unsigned int i, unsigned int j) {
    return (((i * i) - i) >> 1) + j;
} // ltsf0

inline void rltsf0(unsigned int n, unsigned int* i, unsigned int* j) {
    *i = (unsigned int)(0.5 + sqrt(2 * n + 2));
    *j = n - ((*i * *i - *i) >> 1);
} // rltsf0


// lower-triangle matrix representation (1)
//
//    0  1  2  3  4  j
//   -- -- -- -- -- --
// 0|
// 1| 0
// 2| 1  2
// 3| 5  4  3
// 4| 6  7  8  9
// i|         n
//
// i = floor(1/2 + sqrt(2n + 2))
// j0 = n - (i * i - i) / 2
// if (i mod 2 == 1) j = i - j0 - 1 else j = j0

inline void rltsf1(unsigned int n, unsigned int* i, unsigned int* j) {
    *i = (unsigned int)(0.5 + sqrt(2 * n + 2));
    *j = n - ((*i * *i - *i) >> 1);
    if (*i % 2 == 1) *j = *i - *j - 1;
} // rltsf1


// lower-triangle matrix representation (2)
//
//    0  1  2  3  4  j
//   -- -- -- -- -- --
// 0| 0
// 1| 1  2
// 2| 3  4  5
// 3| 6  7  8  9
// i|         n
//
// n = ((i * i) - i) / 2 + i + j

inline unsigned int ltsf2(unsigned int i, unsigned int j) {
    return (((i * i) - i) >> 1) + i + j;
} // ltsf2


// lower-triangle matrix representation (3)
//
//    0  1  2  3  4  j
//   -- -- -- -- -- --
// 0| 0
// 1| 1  2
// 2| 5  4  3
// 3| 6  7  8  9
// i|         n
//
// n = ((i * i) - i) / 2 + i + j

inline unsigned int ltsf3(unsigned int i, unsigned int j) {
    unsigned int b = (((i * i) - i) >> 1) + i;
    return (i % 2) ? b + j : b + i - j;
} // ltsf3


// rectangular matrix representation (0)
//
//    0  1  2  3  4  j     N
//   -- -- -- -- -- -- -- --
// 0| 0  1  2  3  4  5     N-1
// 1| 2N-1           N+1   N
// i|        n
//  |
// M|
//
// i = floor(n / N)
// j0 = n % N
// if (i mod 2 == 1) j = N - j0 else j = j0

inline void rrecsf1(unsigned int n, unsigned int* i, unsigned int* j,
		    unsigned int N) {
    *i = n / N;
    *j = n % N;
    if (*i % 2 == 1) *j = N - *j - 1;
} // rrecsf1


// Z-curve in 2D
// Works only for blocks of power of 2
// This code is not optimized!!!

inline unsigned int zsf(unsigned int z, unsigned int i, unsigned int j) {
    unsigned int r = 0;
    int b = 0;

    bool b0;
    bool b1;

    for (; b < z; ++b) {
	b0 = (i & 1);
	b1 = (j & 1);
	r |= (b1 << (2 * b));
	r |= (b0 << (2 * b + 1));
	i = i >> 1;
	j = j >> 1;
    }

    return r;
} // zsf

inline void rzsf(unsigned int z, unsigned int n, unsigned int* i, unsigned int* j) {
    unsigned int mul = 1;
    int b = 0;
    *i = 0;
    *j = 0;
    for (; b < z; ++b) {
	*j += ((n & 1) * mul);
	n >>= 1;
	*i += ((n & 1) * mul);
	n >>= 1;
	mul <<= 1;
    }
} // rzsf

#endif // SFFUN_H
