/***
 *  $Id$
 **
 *  File: elastic-finalize-seq.cpp
 *  Created: Apr 11, 2013
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
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
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
	std::cout << "Usage: " << ELASTIC_FINALIZE_SHORT << " --input name --output name --map name\n";
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
    AppLog() : argv(), wtime(0), clusters(0), maxc(0), meanc(0), medianc(0) {
	time_t t;
	time(&t);
	date = ctime(&t);
    } // AppLog

    std::string date;
    std::string argv;
    double wtime;
    unsigned int clusters;
    unsigned int minc;
    unsigned int maxc;
    unsigned int meanc;
    unsigned int medianc;

    friend std::ostream& operator<<(std::ostream& os, const AppLog& log) {
	os << "execution date: " << log.date;
	os << "program version: " << ELASTIC_FINALIZE_SHORT << " " << ELASTIC_FINALIZE_VERSION << "\n";
	os << "program options: " << log.argv << "\n";
	os << "walltime used: " << log.wtime << "\n";
	os << "extracted clusters: " << log.clusters << "\n";
	os << "smallest cluster: " << log.minc << "\n";
	os << "largest cluster: " << log.maxc << "\n";
	os << "average cluster: " << log.meanc << "\n";
	os << "median cluster: " << log.medianc << "\n";
	return os;
    } // operator<<
}; // struct AppLog


void welcome() {
    std::cout << ELASTIC_FINALIZE_FULL << " Version " << ELASTIC_FINALIZE_VERSION << "\n";
    std::cout << ELASTIC_FINALIZE_COPYRIGHT << "\n";
    std::cout << "\n";
} // welcome


unsigned int write_plain_cluster(std::ostream& os, const SequenceMap& smap, unsigned int cid,
				 const std::vector<unsigned int>& nodes) {
    unsigned int tot = 0;
    unsigned int n = nodes.size();

    os << ">Cluster " << cid << "\n";

    for (unsigned int i = 0; i < n; ++i) {
	SequenceMap::const_iterator iter = smap.find(nodes[i]);
	if (iter == smap.end()) return 0;

	unsigned int m = iter->name.size();
	tot += m;

	for (unsigned j = 0; j < m; ++j) {
	    os << iter->name[j] << "\n";
	}
    } // for i

    os << std::endl;

    return tot;
} // write_plain_cluster


std::pair<bool, std::string> run(const AppConfig& opt, AppLog& log, Reporter& report) {
    double t0 = get_time();

    std::ifstream fin(opt.input.c_str());
    if (!fin) return std::make_pair(false, "unable to open " + opt.input);

    std::ofstream flog((opt.output + ".eflog").c_str());
    if (!flog) return std::make_pair(false, "unable to create " + opt.output + ".eflog");

    std::ofstream fcls((opt.output + ".clust").c_str());
    if (!flog) return std::make_pair(false, "unable to create " + opt.output + ".clust");

    report << step << "reading map..." << std::endl;

    SequenceMap smap;
    if (smap.read(opt.map) == false) return std::make_pair(false, "unable to read " + opt.map);

    report << step << "formatting clusters, be patient..." << std::endl;

    // read clusters and format output
    jaz::getline_iterator<> gi(fin), end;

    std::vector<std::string> buf;
    std::vector<unsigned int> nodes;

    using namespace boost::accumulators;

    unsigned int nclust = 0;
    accumulator_set<unsigned int, stats<tag::min, tag::max, tag::mean, tag::median> > acc;

    for (; gi != end; ++gi) {
	buf.clear();
	jaz::split(' ', *gi, std::back_inserter(buf));

	unsigned int n = buf.size();
	nodes.resize(n);

	for (unsigned int i = 0; i < n; ++i) {
	    try {
		nodes[i] = boost::lexical_cast<unsigned int>(buf[i]);
	    } catch (boost::bad_lexical_cast& ex) {
		return std::make_pair(false, "incorrect cluster " + boost::lexical_cast<std::string>(nclust));
	    }
	}

	unsigned int tot = write_plain_cluster(fcls, smap, nclust, nodes);
	if (tot == 0) return std::make_pair(false, "error in cluster " + boost::lexical_cast<std::string>(nclust));

	acc(tot);
	nclust++;
    } // for gi

    fcls.close();
    fin.close();

    log.wtime = get_time() - t0;
    log.clusters = nclust;
    log.minc = min(acc);
    log.maxc = max(acc);
    log.meanc = mean(acc);
    log.medianc = median(acc);

    report << info << "extracted " << nclust << " clusters" << std::endl;
    report << info << "smallest cluster: " << log.minc << std::endl;
    report << info << "largest cluster: " << log.maxc << std::endl;
    report << info << "average cluster: " << log.meanc << std::endl;
    report << info << "median cluster: " << log.medianc << std::endl;

    // write final log
    flog << log;
    flog << "config:" << std::endl;
    flog << opt;
    flog.close();

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
