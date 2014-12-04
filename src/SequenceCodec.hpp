/***
 *  $Id$
 **
 *  File: SequenceCodec.hpp
 *  Created: Sep 22, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *          Xiao Yang <eks.yang@gmail.com>
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of CLOSET.
 *  This file is part of ELaSTIC.
 */

#ifndef SEQUENCE_CODEC_HPP
#define SEQUENCE_CODEC_HPP

#include <cstring>
#include <string>


class SequenceCodec {
public:
    SequenceCodec(bool dna = true) : dna_(dna) {
        const char L[4] = { 'A', 'C', 'G', 'T' };

        std::memset(c2s_, '?', 512);
        std::memset(s2c_, '?', 65536);

        for (unsigned int i = 0; i < 16; ++i) {
            unsigned int x = i >> 2;
            unsigned int y = i % 4;
            s2c_[(L[x] << 8) + L[y]] = ('a' + i);
            c2s_[('a' + i) << 1] = L[x];
            c2s_[(('a' + i) << 1) + 1] = L[y];
        }
    }; // SequenceCodec


    const std::string& code(const std::string& s) {
        if (dna_ == true) {
            unsigned int l = s.size();
            unsigned int cl = (l >> 1);

            R_.resize(cl + (l % 2), ' ');

            for (unsigned int i = 0; i < cl; ++i) {
                R_[i] = s2c_[(s[i << 1] << 8) + s[(i << 1) + 1]];
            }

            // if the last letter in odd string is not from [ACGT]
            // we will generate incorrect encoding, in most cases
            // it will result in clear ?? in decoding :-)
            if (l % 2) R_[cl] = s[l - 1];

            return R_;
        }

        return s;
    } // code


    const std::string& decode(const std::string& s) {
        if (dna_ == true) {
            unsigned int l = s.size();
            char c = s[l - 1];

            unsigned int cl = l << 1;
            if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) cl--;

            R_.resize(cl, ' ');

            for (unsigned int i = 0; i < l; ++i) {
                c = s[i];
                R_[i << 1] = c2s_[c << 1];
                R_[(i << 1) + 1] = c2s_[(c << 1) + 1];
            }

            c = s[l - 1];
            if ((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) R_[cl - 1] = c;

            return R_;
        }

        return s;
    } // decode


private:
    bool dna_;

    char c2s_[512];
    char s2c_[65536];

    std::string R_;

}; // SequenceCodec

#endif // SEQUENCE_CODEC_HPP
