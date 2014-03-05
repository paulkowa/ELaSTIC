/***
 *  $Id$
 **
 *  File: tools.hpp
 *  Created: Aug 17, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <sys/time.h>
#include <sys/stat.h>


inline unsigned long int file_size(const char* name) {
    struct stat buf;
    int res = stat(name, &buf);
    return (res == 0 ? static_cast<unsigned long int>(buf.st_size) : 0);
} // file_size

inline double get_time() {
    timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + (0.000001 * t.tv_usec);
} // get_time

#endif // TOOLS_HPP
