/***
 *  $Id$
 **
 *  File: CompressedAlphabet.hpp
 *  Created: Dec 22, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef COMPRESSED_ALPHABET_HPP
#define COMPRESSED_ALPHABET_HPP

#include <cstring>
#include <string>
#include <vector>

#include <jaz/string.hpp>


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

bool any(std::string sigma, char map[]) {
    if (sigma.empty() == true) return false;
    if ((sigma[0] != '[') || (sigma[sigma.size() - 1] != ']')) return false;

    std::vector<std::string> A;

    sigma.erase(sigma.begin());
    sigma.erase(sigma.end() - 1);
    jaz::split(',', sigma, std::back_inserter(A));

    for (unsigned int j = 0; j < A.size(); ++j) {
	std::string& s = A[j];
	unsigned int l = s.size();
	char c = std::toupper(s[0]);
	for (unsigned int i = 0; i < l; ++i) {
	    char c0 = std::tolower(s[i]);
	    char c1 = std::toupper(s[i]);
	    map[c0] = map[c1] = c;
	}
    } // for j

    return true;
} // any

class CompressedAlphabet {
public:
    explicit CompressedAlphabet(const std::string& sigma) {
	std::memset(map_, ' ', 256);
	if (sigma == "A20") A20(map_);
	else if (sigma == "Dayhoff6") Dayhoff6(map_);
	else any(sigma, map_);
    } // CompressedAlphabet

    bool test() const {
	if ((map_['A'] == ' ') || (map_['a'] == ' ')) return false;
	if ((map_['C'] == ' ') || (map_['c'] == ' ')) return false;
	if ((map_['D'] == ' ') || (map_['d'] == ' ')) return false;
	if ((map_['E'] == ' ') || (map_['e'] == ' ')) return false;
	if ((map_['F'] == ' ') || (map_['f'] == ' ')) return false;
	if ((map_['G'] == ' ') || (map_['g'] == ' ')) return false;
	if ((map_['H'] == ' ') || (map_['h'] == ' ')) return false;
	if ((map_['I'] == ' ') || (map_['i'] == ' ')) return false;
	if ((map_['K'] == ' ') || (map_['k'] == ' ')) return false;
	if ((map_['L'] == ' ') || (map_['l'] == ' ')) return false;
	if ((map_['M'] == ' ') || (map_['m'] == ' ')) return false;
	if ((map_['N'] == ' ') || (map_['n'] == ' ')) return false;
	if ((map_['P'] == ' ') || (map_['p'] == ' ')) return false;
	if ((map_['Q'] == ' ') || (map_['q'] == ' ')) return false;
	if ((map_['R'] == ' ') || (map_['r'] == ' ')) return false;
	if ((map_['S'] == ' ') || (map_['s'] == ' ')) return false;
	if ((map_['T'] == ' ') || (map_['t'] == ' ')) return false;
	if ((map_['V'] == ' ') || (map_['v'] == ' ')) return false;
	if ((map_['W'] == ' ') || (map_['w'] == ' ')) return false;
	if ((map_['Y'] == ' ') || (map_['y'] == ' ')) return false;
	return true;
    } // test

    std::string operator()(std::string s) const {
	unsigned int l = s.size();
	for (unsigned int i = 0; i < l; ++i) s[i] = map_[s[i]];
	return s;
    } // operator()

private:
    char map_[256];

}; // CompressedAlphabet

#endif // COMPRESSED_ALPHABET_HPP
