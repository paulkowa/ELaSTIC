/***
 *  $Id$
 **
 *  File: elastic-sketch-omp.cpp
 *  Created: Nov 17, 2014
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2014 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#include <iostream>
#include <map>
#include <new>
#include <string>
#include <sstream>

#include "config.hpp"
#include "iomanip.hpp"

#include "sketch/config_log.hpp"

#include <jaz/parameters.hpp>
#include <jaz/hash.hpp>


void welcome() {
    std::cout << "ELASTIC_SKETCH_FULL" << " Version " << "ELASTIC_SKETCH_VERSION" << "\n";
    std::cout << "ELASTIC_SKETCH_COPYRIGHT" << "\n";
    std::cout << std::endl;
} // welcome


void run(const AppConfig& opt, AppLog& log, Reporter& report) {
    bool res = true;
    std::string err = "";

} // run


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    welcome();

    // get parameters
    std::map<std::string, std::string> conf;

    if (argc == 1) {
	AppConfig::usage();
	return 0;
    }

    bool res = false;
    int pos = 0;

    boost::tie(res, pos) = jaz::parse_argv(argc, argv, conf);

    if (res == false) {
	if (pos == -1) {
	    AppConfig::usage();
	    std::cout << error << "incorrect command line arguments\n";
	    return -1;
	} else {
	    std::cout << error << "incorrect command line argument " << argv[pos] << "\n";
	    return -1;
	}
    } // if res

    // set reporter
    Reporter report(std::cout, std::cout);

    // create config and log
    AppConfig opt;
    std::string err = "";

    boost::tie(res, err) = opt.set(conf);

    if (res == false) {
	report << error << err << "\n";
	return -1;
    }

    // TODO: update for OpenMP
    AppLog log;
    log.argv = jaz::join(' ', argv + 1, argv + argc);
    log.cpus = 0;

    // and here we go
    run(opt, log, report);

    return 0;
} // main
