/***
 *  $Id$
 **
 *  File: write_output.hpp
 *  Created: Mar 21, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef WRITE_OUTPUT_HPP
#define WRITE_OUTPUT_HPP

#include <iterator>
#include <string>
#include <sstream>

#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"
#include "tools.hpp"

#include <mpix2/write_cbuffer.hpp>


inline std::pair<bool, std::string> write_output(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
                                                 const std::vector<read_pair>& edges, unsigned long int& fs) {
    report << info << "writing output..." << std::endl;

    const int MPI_ABRT_SIG = 13;

    /*
      boost::format fmt("%05d");
      fmt % rank;

      std::ofstream of((opt.output + ".sim." + fmt.str()).c_str());
      std::copy(edges.begin(), edges.end(), std::ostream_iterator<read_pair>(of, "\n"));
      of.close();
    */

    std::string name;
    std::ostringstream os;

    if (opt.factor == true) {
        name = opt.output + ".tim.00000";
        unsigned int m = edges.size();
        for (unsigned int i = 0; i < m; ++i) {
            write_read_pair(os, edges[i]);
            os << "\n";
        }
    } else {
        name = opt.output + ".sim.00000";
        std::copy(edges.begin(), edges.end(), std::ostream_iterator<read_pair>(os, "\n"));
    }

    if (mpix::write_cbuffer(name, os.str().c_str(), os.str().size(), comm) == false) {
        report.critical << error << ("unable to create " + name) << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }

    fs = file_size(name.c_str());

    return std::make_pair(true, "");
} // write_output

#endif // WRITE_OUTPUT_HPP
