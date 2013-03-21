/***
 *  $Id$
 **
 *  File: write_output.hpp
 *  Created: Mar 21, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
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

#include <mpix/write_cbuffer.hpp>


inline std::pair<bool, std::string> write_output(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm,
						 const std::vector<read_pair>& edges) {
    report << info << "writing output..." << std::endl;

    const int MPI_ABRT_SIG = 13;

    /*
      boost::format fmt("%05d");
      fmt % rank;

      std::ofstream of((opt.output + "." + fmt.str()).c_str());
      std::copy(edges.begin(), edges.end(), std::ostream_iterator<read_pair>(of, "\n"));
      of.close();
    */

    std::ostringstream os;
    std::copy(edges.begin(), edges.end(), std::ostream_iterator<read_pair>(os, "\n"));

    if (mpix::write_cbuffer((opt.output + ".00000"), os.str().c_str(), os.str().size(), comm) == false) {
	report.critical << error << ("unable to create " + opt.output + ".00000") << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }

    return std::make_pair(true, "");
} // write_output

#endif // WRITE_OUTPUT_HPP
