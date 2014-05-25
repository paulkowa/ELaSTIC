/***
 *  $Id$
 **
 *  File: iomanip.hpp
 *  Created: May 20, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef IOMANIP_HPP
#define IOMANIP_HPP

#include <iostream>


struct null_ostream : public std::ostream {
    null_ostream() : std::ios(0), std::ostream(0) { }
}; // struct null_ostream

extern null_ostream nout;


struct Reporter {
    Reporter() : stream(std::cout), normal(std::cout), critical(std::cerr) { }
    Reporter(std::ostream& out) : stream(out), normal(out), critical(out) { }
    Reporter(std::ostream& out, std::ostream& err) : stream(out), normal(out), critical(err) { }
    Reporter(std::ostream& str, std::ostream& out, std::ostream& err) : stream(str), normal(out), critical(err) { }

    std::ostream& stream;
    std::ostream& normal;
    std::ostream& critical;
}; // struct Reporter

inline std::ostream& operator<<(Reporter& report, std::ostream& (*f)(std::ostream&)) {
    return report.normal << f;
} // operator<<

template <typename T> inline std::ostream& operator<<(Reporter& report, const T& t) {
    return report.normal << t;
} // operator<<


inline std::ostream& warning(std::ostream& os) {
    os << "~ warning: ";
    return os;
} // warning

inline std::ostream& error(std::ostream& os) {
    os << "! error: ";
    return os;
} // error

inline std::ostream& option(std::ostream& os) {
    os << "";
    return os;
} // option

inline std::ostream& info(std::ostream& os) {
    os << "";
    return os;
} // info

inline std::ostream& step(std::ostream& os) {
    os << "";
    return os;
} // info

inline std::ostream& clock(std::ostream& os) {
    os << "";
    return os;
} // clock

inline std::ostream& debug(std::ostream& os) {
    os << "...";
    return os;
} // debug

class timer {
public:
    explicit timer(double t) : t_(t) { }

private:
    double t_;

    friend std::ostream& operator<<(std::ostream& os, const timer& t) {
	unsigned int tt = static_cast<unsigned int>(t.t_);
	unsigned int ht = tt / 3600;
	unsigned int mt = (tt % 3600) / 60;
	unsigned int st = (tt % 3600) % 60;
	os << t.t_ << "s [" << ht << "h" << mt << "m" << st << "s]";
	return os;
    } // operator <<

}; // class timer

#endif // IOMANIP_HPP
