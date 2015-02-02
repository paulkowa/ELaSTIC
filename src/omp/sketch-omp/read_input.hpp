/***
 *  $Id$
 **
 *  File: read_input.hpp
 *  Created: Dec 16, 2014
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2014 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef READ_INPUT_HPP
#define READ_INPUT_HPP

#include <numeric>
#include <arpa/inet.h>

#include "sketch/Sequence.hpp"

#include "SequenceCodec.hpp"
#include "tools.hpp"


inline std::pair<bool, std::string> read_input(const AppConfig& opt, AppLog& log, Reporter& report,
                                               SequenceList& SL) {
    report << step << "reading input sequences..." << std::endl;

    std::vector<Sequence>& seqs = SL.seqs;

    // read index
    unsigned long int fsz = file_size((opt.input + ".eidx").c_str());
    if (fsz == 0) return std::make_pair(false, "unable to open " + opt.input + ".eidx");

    typedef uint16_t index_type;
    std::vector<index_type> index;

    unsigned int n = fsz / sizeof(index_type);
    SL.N = n;

    index.resize(n);
    seqs.resize(n);

    std::ifstream f((opt.input + ".eidx").c_str(), std::ifstream::binary);
    f.read(reinterpret_cast<char*>(&index[0]), fsz);
    f.close();

    std::transform(index.begin(), index.end(), index.begin(), ntohs);

    // read data
    fsz = file_size((opt.input + ".eseq").c_str());
    if (fsz == 0) return std::make_pair(false, "unable to open " + opt.input + ".eseq");

    std::vector<char> data(fsz);

    f.open((opt.input + ".eseq").c_str(), std::ifstream::binary);
    f.read(&data[0], fsz);
    f.close();

    // decode data
    SequenceCodec sc(opt.is_dna);

    unsigned int pos = index[0];

    seqs[0].id = ntohl(*reinterpret_cast<uint32_t*>(&data[0]));
    seqs[0].s = sc.decode(std::string(data.begin() + 4, data.begin() + pos));

    if ((seqs[0].s.size() < opt.kmer) || (seqs[0].s.find('?') != std::string::npos)) {
        return std::make_pair(false, "invalid input sequence, change kmer?");
    }

    for (unsigned int i = 1; i < n; ++i) {
        seqs[i].id = ntohl(*reinterpret_cast<uint32_t*>(&data[pos]));
        seqs[i].s = sc.decode(std::string(data.begin() + pos + 4, data.begin() + pos + index[i]));

        if ((seqs[i].s.size() < opt.kmer) || (seqs[i].s.find('?') != std::string::npos)) {
            return std::make_pair(false, "invalid input sequence, change kmer?");
        }

        pos += index[i];
    }

    // update log
    log.input = n;

    report << info << "found " << n << " sequences" << std::endl;

    return std::make_pair(true, "");
} // read_input

#endif // READ_INPUT_HPP
