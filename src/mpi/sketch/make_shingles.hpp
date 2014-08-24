/***
 *  $Id$
 **
 *  File: make_shingles.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef MAKE_SHINGLES_HPP
#define MAKE_SHINGLES_HPP

#include <string>
#include <vector>

#include "CompressedAlphabet.hpp"
#include "SequenceCodec.hpp"

#include "SequenceDB.hpp"
#include "config_log.hpp"


template <typename Hash>
std::pair<bool, std::string> make_shingles(const AppConfig& opt, AppLog& log, Reporter& report,
					   const SequenceList& SL, const Hash& hash,
					   std::vector<shingle_list_type>& shingles) {
    report << step << "creating shingles..." << std::endl;

    unsigned short int k = opt.kmer;
    unsigned int n = SL.seqs.size();

    CompressedAlphabet ca(opt.sigma);
    if (opt.is_dna == true) ca = CompressedAlphabet("A20");

    SequenceCodec sc(opt.is_dna);
    shingles.resize(n);

    for (unsigned int i = 0; i < n; ++i) {
	unsigned int l = SL.seqs[i].s.size();
	std::string s = ca(SL.seqs[i].s);

	shingles[i].resize(l - k + 1);
	for (unsigned int j = 0; j < l - k + 1; ++j) shingles[i][j] = hash(sc.code(s.substr(j, k)));
    } // for i

    return std::make_pair(true, "");
} // make_shingles

#endif // MAKE_SHINGLES_HPP
