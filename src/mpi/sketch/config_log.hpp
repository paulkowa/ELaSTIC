/***
 *  $Id$
 **
 *  File: config_log.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef CONFIG_LOG_HPP
#define CONFIG_LOG_HPP

#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <jaz/parameters.hpp>


struct AppConfig {
    AppConfig() {
	input = "";
	output = "";
	method = 0;
	level = 75;
	kmer = 15;
	mod = 25;
	iter = -1;
	cmax = 5000;
	jmin = 50;
	rma = 0;
    } // AppConfig

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

	// check config
	Container ext_conf;

	if (jaz::check_option(conf, "config", val) == true) {
	    bool res = true;
	    int line = 0;

	    boost::tie(res, line) = jaz::read_config(val, ext_conf);

	    if (res == false) {
		return std::make_pair(false, "incorrect config file, line " + boost::lexical_cast<std::string>(line));
	    }

	    // update config
	    typename Container::const_iterator iter(conf.begin());
	    for (; iter != conf.end(); ++iter) {
		typename Container::iterator pos;
		boost::tie(pos, res) = ext_conf.insert(*iter);
		if (res == false) pos->second = iter->second;
	    } // for iter
	} else ext_conf = conf;

	// from here complete config is in ext_conf

	if (jaz::check_option(ext_conf, "method", val) == true) {
	    method = boost::lexical_cast<unsigned short int>(val);
	    if (method > 1) {
		return std::make_pair(false, "incorrect method");
	    }
	}

	if (jaz::check_option(ext_conf, "level", val) == true) {
	    level = boost::lexical_cast<short int>(val);
	    if ((level < 10) || (level > 100)) {
		return std::make_pair(false, "incorrect level");
	    }
	}

	if (jaz::check_option(ext_conf, "kmer", val) == true) {
	    kmer = boost::lexical_cast<unsigned int>(val);
	    if ((kmer < 5) || (kmer > 31)) {
		return std::make_pair(false, "incorrect kmer");
	    }
	}

	if (jaz::check_option(ext_conf, "mod", val) == true) {
	    mod = boost::lexical_cast<short int>(val);
	    if ((mod < 1) || (mod > 100)) {
		return std::make_pair(false, "incorrect mod");
	    }
	}

	if (jaz::check_option(ext_conf, "iter", val) == true) {
	    iter = boost::lexical_cast<short int>(val);
	    if ((iter < 1) || (mod < iter)) {
		return std::make_pair(false, "incorrect iter");
	    }
	}

	if (iter == -1) iter = mod;

	if (jaz::check_option(ext_conf, "cmax", val) == true) {
	    cmax = boost::lexical_cast<short int>(val);
	    if (cmax < 100) {
		return std::make_pair(false, "incorrect cmax");
	    }
	}

	if (jaz::check_option(ext_conf, "jmin", val) == true) {
	    jmin = boost::lexical_cast<short int>(val);
	    if ((jmin < 10) || (jmin > 100)) {
		return std::make_pair(false, "incorrect jmin");
	    }
	}

	// temporal options
	if (jaz::check_option(ext_conf, "rma", val) == true) {
	    rma = boost::lexical_cast<short int>(val);
	}

	return std::make_pair(true, "");
    } // set

    std::string input;
    std::string output;
    unsigned short int method;
    unsigned short int level;
    unsigned short int kmer;
    short int mod;
    short int iter;
    unsigned int cmax;
    unsigned short int jmin;
    short int rma;

}; // struct AppConfig


struct AppLog {
    AppLog() : argv(), cpus(0), wtime(0), input(0), cedges(0), vedges(0), through(0) {
	time_t t;
	time(&t);
	date = ctime(&t);
    } // AppLog

    std::string date;
    std::string argv;
    unsigned int cpus;
    double wtime;
    unsigned int input;
    unsigned int cedges;
    unsigned int vedges;
    double through;

    friend std::ostream& operator<<(std::ostream& os, const AppLog& log) {
	os << "execution date: " << log.date;
	os << "program version: " << ELASTIC_SKETCH_SHORT << " " << ELASTIC_SKETCH_VERSION << "\n";
	os << "program options: " << log.argv << "\n";
	os << "processors used: " << log.cpus << "\n";
	os << "walltime used: " << log.wtime << "\n";
	os << "input sequences: " << log.input << "\n";
	os << "candidate edges: " << log.cedges << "\n";
	os << "validated edges: " << log.vedges << "\n";
	os << "edge throughput: " << log.through << "\n";
	return os;
    } // operator<<
}; // struct AppLog

#endif // CONFIG_LOG_HPP
