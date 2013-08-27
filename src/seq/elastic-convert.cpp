/***
 *  $Id$
 **
 *  File: elastic-convert.cpp
 *  Created: Aug 19, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#include <fstream>
#include <iostream>
#include <map>

#include <jaz/iterator.hpp>
#include <jaz/parameters.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "SequenceMap.hpp"
#include "config.hpp"
#include "iomanip.hpp"
#include "tools.hpp"


struct AppConfig {
    AppConfig() {
	input = "";
	output = "";
	map = "";
    } // AppConfig

    static void usage() {
	std::cout << "Usage: " << ELASTIC_CONVERT_SHORT << " --input name --output name --map name\n";
	std::cout << "\n";
	std::cout << "Options:\n";
	std::cout << "  --input name          read input from this file\n";
	std::cout << "  --output name         write output to files with this prefix\n";
	std::cout << "  --map name            read sequence map from this file\n";
	std::cout << "\n";
    } // usage

    template <typename Container>
    std::pair<bool, std::string> set(const Container& conf) {
	std::string val;

	// check major options
	if (jaz::check_option(conf, "input", val) == false) {
	    return std::make_pair(false, "missing input parameter");
	}

	input = val;

	if (jaz::check_option(conf, "output", val) == false) {
	    return std::make_pair(false, "missing output parameter");
	}

	output = val;

	if (jaz::check_option(conf, "map", val) == false) {
	    return std::make_pair(false, "missing map parameter");
	}

	map = val;

	return std::make_pair(true, "");
    } // set

    std::string input;
    std::string output;
    std::string map;

    friend std::ostream& operator<<(std::ostream& os, const AppConfig& opt) {
	os << "input = " << opt.input << "\n";
	os << "output = " << opt.output << "\n";
	os << "map = " << opt.map << "\n";
	return os;
    } // operator<<

}; // struct AppConfig


struct AppLog {
    AppLog() : argv(), wtime(0) {
	time_t t;
	time(&t);
	date = ctime(&t);
    } // AppLog

    std::string date;
    std::string argv;
    double wtime;

    friend std::ostream& operator<<(std::ostream& os, const AppLog& log) {
	os << "execution date: " << log.date;
	os << "program version: " << ELASTIC_CONVERT_SHORT << " " << ELASTIC_CONVERT_VERSION << "\n";
	os << "program options: " << log.argv << "\n";
	os << "walltime used: " << log.wtime << "\n";
	return os;
    } // operator<<
}; // struct AppLog


void welcome() {
    std::cout << ELASTIC_CONVERT_FULL << " Version " << ELASTIC_CONVERT_VERSION << "\n";
    std::cout << ELASTIC_CONVERT_COPYRIGHT << "\n";
    std::cout << "\n";
} // welcome


std::pair<bool, std::string> run(const AppConfig& opt, AppLog& log, Reporter& report) {
    double t0 = get_time();

    std::ifstream fin(opt.input.c_str());
    if (!fin) return std::make_pair(false, "unable to open " + opt.input);

    std::ofstream fout((opt.output + ".tsv").c_str());
    if (!fout) return std::make_pair(false, "unable to create " + opt.output + ".tsv");

    report << step << "reading map..." << std::endl;

    SequenceMap smap;
    if (smap.read(opt.map) == false) return std::make_pair(false, "unable to read " + opt.map);

    report << step << "converting graph, be patient..." << std::endl;

    jaz::getline_iterator<> iter(fin), end;

    for (; iter != end; ++iter) {
	std::istringstream is(*iter);

	unsigned int s;
	is >> s;
	if (!is) return std::make_pair(false, "unable to read " + opt.input);

	unsigned int t;
	is >> t;
	if (!is) return std::make_pair(false, "unable to read " + opt.input);

	std::string ln;
	std::getline(is, ln);
	std::replace(ln.begin(), ln.end(), ' ', '\t');

	SequenceMap::const_iterator siter0 = smap.find(s);
	SequenceMap::const_iterator siter1 = smap.find(t);

	if ((siter0 == smap.end()) || (siter1 == smap.end())) return std::make_pair(false, "incorrect index");

	fout << siter0->name[0] << "\t" << siter1->name[0] << ln << std::endl;
    } // for iter

    fout.close();
    log.wtime = get_time() - t0;

    return std::make_pair(true, "");
} // run


int main(int argc, char* argv[]) {
    welcome();

    // get parameters
    std::map<std::string, std::string> conf;

    if (argc == 1) {
	AppConfig::usage();
	return 0;
    }

    bool res = false;
    int pos = 0;
    std::string err = "";

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
    }

    // create config and log
    AppConfig opt;
    boost::tie(res, err) = opt.set(conf);

    if (res == false) {
	std::cout << error << err << "\n";
	return -1;
    }

    AppLog log;
    log.argv = jaz::join(' ', argv + 1, argv + argc);

    Reporter report(std::cout, std::cout);

    boost::tie(res, err) = run(opt, log, report);

    if (res == false) {
	std::cout << error << err << "\n";
	return -1;
    }

    std::cout << "time: " << log.wtime << std::endl;
    std::cout << "done!" << std::endl;
    std::cout << std::endl;

    return 0;
} // main