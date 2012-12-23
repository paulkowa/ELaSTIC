/***
 *  $Id$
 **
 *  File: CompressedAlphabet.hpp
 *  Created: Dec 22, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef COMPRESSEDALPHABET_HPP
#define COMPRESSEDALPHABET_HPP

#include <cstring>
#include <string>


void A20(char map[]) {
    map['A'] = map['a'] = 'A';
    map['C'] = map['c'] = 'C';
    map['D'] = map['d'] = 'D';
    map['E'] = map['e'] = 'E';
    map['F'] = map['f'] = 'F';
    map['G'] = map['g'] = 'G';
    map['H'] = map['h'] = 'H';
    map['I'] = map['i'] = 'I';
    map['K'] = map['k'] = 'K';
    map['L'] = map['l'] = 'L';
    map['M'] = map['m'] = 'M';
    map['N'] = map['n'] = 'N';
    map['P'] = map['p'] = 'P';
    map['Q'] = map['q'] = 'Q';
    map['R'] = map['r'] = 'R';
    map['S'] = map['s'] = 'S';
    map['T'] = map['t'] = 'T';
    map['V'] = map['v'] = 'V';
    map['W'] = map['w'] = 'W';
    map['Y'] = map['y'] = 'Y';
} // A20

void Dayhoff6(char map[]) {
    map['A'] = map['a'] = 'A';
    map['C'] = map['c'] = 'C';
    map['D'] = map['d'] = 'D';
    map['E'] = map['e'] = 'D';
    map['F'] = map['f'] = 'F';
    map['G'] = map['g'] = 'A';
    map['H'] = map['h'] = 'H';
    map['I'] = map['i'] = 'I';
    map['K'] = map['k'] = 'H';
    map['L'] = map['l'] = 'I';
    map['M'] = map['m'] = 'I';
    map['N'] = map['n'] = 'D';
    map['P'] = map['p'] = 'A';
    map['Q'] = map['q'] = 'D';
    map['R'] = map['r'] = 'H';
    map['S'] = map['s'] = 'A';
    map['T'] = map['t'] = 'A';
    map['V'] = map['v'] = 'I';
    map['W'] = map['w'] = 'F';
    map['Y'] = map['y'] = 'F';
} // Dayhoff6


class CompressedAlphabet {
public:
    explicit CompressedAlphabet(const std::string& sigma) {
	std::memset(map_, ' ', 256);
	if (sigma == "A20") A20(map_);
	else if (sigma == "Dayhoff6") Dayhoff6(map_);
    } // CompressedAlphabet

    std::string operator()(std::string s) const {
	unsigned int l = s.size();
	for (unsigned int i = 0; i < l; ++i) s[i] = map_[s[i]];
	return s;
    } // operator()

private:
    char map_[256];

}; // CompressedAlphabet

#endif // COMPRESSEDALPHABET_HPP
