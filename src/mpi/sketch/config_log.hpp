/***
 *  $Id$
 **
 *  File: config_log.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef CONFIG_LOG_HPP
#define CONFIG_LOG_HPP

#include <iomanip>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <jaz/parameters.hpp>

#include "CompressedAlphabet.hpp"
#include "create_smatrix.hpp"

#include "config.hpp"


struct AppConfig {
    AppConfig() {
	input = "";
	output = "";
	is_dna = true;
	sigma = "A20";
	compress = 0;
	method = 0;
	kmer = 16;
	gaps = "[1,-2,-10,-1]";
	level = 75;
	factor = false;
	mod = 25;
	iter = 7;
	cmax = 10000;
	jmin = 50;
	wsq = true;
	dbg = -1;
	eps = 0;

	params.push_back("input");
	params.push_back("output");
	params.push_back("config");
	params.push_back("type");
	params.push_back("sigma");
	params.push_back("compress");
	params.push_back("method");
	params.push_back("kmer");
	params.push_back("gaps");
	params.push_back("level");
	params.push_back("factor");
	params.push_back("modulo");
	params.push_back("iterate");
	params.push_back("cmax");
	params.push_back("jmin");
	params.push_back("steal");
	params.push_back("debug");
	params.push_back("eps");
    } // AppConfig

    static void usage() {
	std::cout << "Usage: " << ELASTIC_SKETCH_SHORT << " --input name --output name [options...]\n";
	std::cout << "\n";
	std::cout << "Options:\n";
	std::cout << "  --input name       read input from files with this prefix\n";
	std::cout << "  --output name      write output to files with this prefix\n";
	std::cout << "  --config name      read configuration from this file\n";
	std::cout << "  --type {nt|aa}     set input sequence type (default nt)\n";
	std::cout << "  --sigma type       use this compressed amino acid alphabet (default A20)\n";
	std::cout << "  --compress {0|1}   use compressed alphabet during validation (default 0)\n";
	std::cout << "  --method type      use this method for validation (default 0)\n";
	std::cout << "  --kmer size        use kmers of this size (default 16)\n";
	std::cout << "  --gaps type        use these alignment parameters (default [1,-2,-10,-1])\n";
	std::cout << "  --level size       use this threshold during validation (default 75)\n";
	std::cout << "  --factor {0|1}     output intermediate values of similarity score (default 0)\n";
	std::cout << "  --modulo size      use this mod value in sketching (default 25)\n";
	std::cout << "  --iterate size     limit the number of sketching iterations to this (default 7)\n";
	std::cout << "  --cmax size        use this limit to mark frequent kmers (default 10000)\n";
	std::cout << "  --jmin size        use this limit to extract candidate pairs (default 50)\n";
	std::cout << "  --steal {0|1}      enable work stealing during validation (default 1)\n";
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

	try {

	    if (jaz::check_option(ext_conf, "method", val) == true) {
		method = boost::lexical_cast<unsigned short int>(val);
		if (method > 5) {
		    return std::make_pair(false, "incorrect method");
		}
	    }

	    if (jaz::check_option(ext_conf, "type", val) == true) {
		if ((val != "nt") && (val != "aa")) {
		    return std::make_pair(false, "incorrect type");
		}
		is_dna = (val == "nt");
	    }

	    if (jaz::check_option(ext_conf, "sigma", val) == true) {
		if ((val == "A20") || (val == "Dayhoff6")) sigma = val;
		else {
		    CompressedAlphabet ca(val);
		    if (ca.test() == false) return std::make_pair(false, "unknown compressed alphabet");
		    sigma = val;
		}
	    }

	    if (jaz::check_option(ext_conf, "compress", val) == true) {
		compress = boost::lexical_cast<bool>(val);
	    }

	    if (jaz::check_option(ext_conf, "kmer", val) == true) {
		kmer = boost::lexical_cast<unsigned int>(val);
		if ((kmer < 3) || (kmer > 31)) {
		    return std::make_pair(false, "incorrect kmer");
		}
	    }

	    if (jaz::check_option(ext_conf, "gaps", val) == true) {
		gaps = val;

		int g, h;
		bio::scoring_matrix sm;

		if (create_smatrix(gaps, is_dna, sm, g, h) == false) {
		    return std::make_pair(false, "incorrect gaps");
		}
	    }

	    if (jaz::check_option(ext_conf, "level", val) == true) {
		level = boost::lexical_cast<short int>(val);
		if (level > 100) {
		    return std::make_pair(false, "incorrect level");
		}
	    }

	    if (jaz::check_option(ext_conf, "factor", val) == true) {
		factor = boost::lexical_cast<bool>(val);
	    }

	    if (jaz::check_option(ext_conf, "modulo", val) == true) {
		mod = boost::lexical_cast<short int>(val);
		if ((mod < 1) || (mod > 100)) {
		    return std::make_pair(false, "incorrect modulo");
		}
	    }

	    if (jaz::check_option(ext_conf, "iterate", val) == true) {
		iter = boost::lexical_cast<short int>(val);
		if ((iter < 1) || (mod < iter)) {
		    return std::make_pair(false, "incorrect iterate");
		}
	    }

	    if (jaz::check_option(ext_conf, "cmax", val) == true) {
		cmax = boost::lexical_cast<unsigned int>(val);
		if (cmax < 1) {
		    return std::make_pair(false, "incorrect cmax");
		}
	    }

	    if (jaz::check_option(ext_conf, "jmin", val) == true) {
		jmin = boost::lexical_cast<short int>(val);
		if ((jmin < 10) || (jmin > 100)) {
		    return std::make_pair(false, "incorrect jmin");
		}
	    }

	    if (jaz::check_option(ext_conf, "steal", val) == true) {
		wsq = boost::lexical_cast<bool>(val);
	    }

	    if (jaz::check_option(ext_conf, "debug", val) == true) {
		dbg = boost::lexical_cast<int>(val);
		if (dbg < 0) {
		    return std::make_pair(false, "incorrect debug");
		}
	    }

	    if (jaz::check_option(ext_conf, "eps", val) == true) {
		eps = boost::lexical_cast<unsigned int>(val);
		if (eps > cmax) {
		    return std::make_pair(false, "incorrect eps");
		}
	    }
	} catch (boost::bad_lexical_cast& ex) {
	    return std::make_pair(false, "incorrect type of argument(s)");
	}

	for (unsigned int i = 0; i < params.size(); ++i) {
	    typename Container::iterator iter(ext_conf.find(params[i]));
	    if (iter != ext_conf.end()) ext_conf.erase(iter);
	}

	if (ext_conf.empty() == false) {
	    return std::make_pair(false, "unknown parameter " + ext_conf.begin()->first);
	}

	return std::make_pair(true, "");
    } // set

    std::string input;
    std::string output;
    bool is_dna;
    std::string sigma;
    bool compress;
    unsigned short int method;
    unsigned short int kmer;
    std::string gaps;
    unsigned short int level;
    bool factor;
    short int mod;
    short int iter;
    unsigned int cmax;
    unsigned short int jmin;
    bool wsq;
    int dbg; // internal: print debug information
    unsigned int eps; // internal: partition size in sketching balance

    std::vector<std::string> params;

    friend std::ostream& operator<<(std::ostream& os, const AppConfig& opt) {
	os << "input = " << opt.input << "\n";
	os << "output = " << opt.output << "\n";
	os << "dna = " << opt.is_dna << "\n";
	os << "sigma = " << opt.sigma << "\n";
	os << "compress = " << opt.compress << "\n";
	os << "method = " << opt.method << "\n";
	os << "kmer = " << opt.kmer << "\n";
	os << "gaps = " << opt.gaps << "\n";
	os << "level = " << opt.level << "\n";
	os << "factor = " << opt.factor << "\n";
	os << "modulo = " << opt.mod << "\n";
	os << "iterate = " << opt.iter << "\n";
	os << "cmax = " << opt.cmax << "\n";
	os << "jmin = " << opt.jmin << "\n";
	os << "steal = " << opt.wsq << "\n";
	return os;
    } // operator<<

}; // struct AppConfig


struct AppLog {
    AppLog() : argv(), cpus(0), wtime(0), input(0), cedges(0), vedges(0) {
	time_t t;
	time(&t);
	date = ctime(&t);
    } // AppLog

    std::string date;
    std::string argv;
    unsigned int cpus;
    double wtime;
    unsigned int input;
    unsigned long long int cedges;
    unsigned long long int vedges;

    friend std::ostream& operator<<(std::ostream& os, const AppLog& log) {
	os << "execution date: " << log.date;
	os << "program version: " << ELASTIC_SKETCH_SHORT << " " << ELASTIC_SKETCH_VERSION << "\n";
	os << "program options: " << log.argv << "\n";
	os << "processors used: " << log.cpus << "\n";
	os << "walltime used: " << log.wtime << "\n";
	os << "input sequences: " << log.input << "\n";
	os << "candidate edges: " << log.cedges << "\n";
	os << "validated edges: " << log.vedges << "\n";
	return os;
    } // operator<<
}; // struct AppLog

#endif // CONFIG_LOG_HPP
