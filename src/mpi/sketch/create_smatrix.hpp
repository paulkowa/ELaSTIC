/***
 *  $Id$
 **
 *  File: create_smatrix.hpp
 *  Created: Jun 20, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef CREATE_SMATRIX_HPP
#define CREATE_SMATRIX_HPP

#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <bio/sequence_compare.hpp>
#include <jaz/string.hpp>


inline bool create_smatrix(std::string val, bool is_dna, bio::scoring_matrix& sm, int& g, int& h) {
    val.erase(val.begin());
    val.erase(val.end() - 1);

    std::vector<std::string> agap;
    jaz::split(',', val, std::back_inserter(agap));

    try {
        if (agap.size() == 3) {
            if (bio::read_file_sm(agap[0], sm) == false) return false;
            g = boost::lexical_cast<int>(agap[1]);
            h = boost::lexical_cast<int>(agap[2]);
        } else if (agap.size() == 4) {
            int gaps[4];
            for (unsigned int i = 0; i < 4; ++i) {
                gaps[i] = boost::lexical_cast<int>(agap[i]);
            }
            if (is_dna == true) sm = bio::make_dna_sm(gaps[0], gaps[1]);
            else sm = bio::make_dummy_sm(gaps[0], gaps[1]);
            g = gaps[2];
            h = gaps[3];
        } else return false;
    } catch (boost::bad_lexical_cast& ex) {
        return false;
    }

    return true;
} // create_smatrix

#endif // CREATE_SMATRIX_HPP
