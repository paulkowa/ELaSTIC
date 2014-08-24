/***
 *  $Id$
 **
 *  File: elastic-cluster-seq.cpp
 *  Created: Apr 05, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#include <fstream>
#include <iostream>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <jaz/parameters.hpp>
#include <jaz/boost/files.hpp>
#include <jaz/science/disjoint_set.hpp>

#include "config.hpp"
#include "iomanip.hpp"
#include "tools.hpp"


struct AppConfig {
    AppConfig() {
	input = "";
	output = "";
	nodes = 0;
	level = 0;
    } // AppConfig

    static void usage() {
	std::cout << "Usage: " << ELASTIC_CLUSTER_SHORT << " --input name --output name --nodes size [options...]\n";
	std::cout << "\n";
	std::cout << "Options:\n";
	std::cout << "  --input name          read input from this file/directory\n";
	std::cout << "  --output name         write output to files with this prefix\n";
	std::cout << "  --nodes size          assume that many input nodes\n";
	std::cout << "  --level size          use this threshold to filter edges (default 0)\n";
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

	try {
	    if (jaz::check_option(conf, "nodes", val) == true) {
		nodes = boost::lexical_cast<unsigned int>(val);
	    }

	    if (jaz::check_option(conf, "level", val) == true) {
		level = boost::lexical_cast<unsigned int>(val);
	    }
	} catch (boost::bad_lexical_cast& ex) {
	    return std::make_pair(false, "incorrect argument(s)");
	}

	if (nodes < 2) {
	    return std::make_pair(false, "incorrect number of nodes");
	}

	return std::make_pair(true, "");
    } // set

    std::string input;
    std::string output;
    unsigned int nodes;
    unsigned int level;

    friend std::ostream& operator<<(std::ostream& os, const AppConfig& opt) {
	os << "input = " << opt.input << "\n";
	os << "output = " << opt.output << "\n";
	os << "nodes = " << opt.nodes << "\n";
	os << "level = " << opt.level << "\n";
	return os;
    } // operator<<

}; // struct AppConfig


struct AppLog {
    AppLog() : argv(), wtime(0), input(0), extracted(0) {
	time_t t;
	time(&t);
	date = ctime(&t);
    } // AppLog

    std::string date;
    std::string argv;
    double wtime;
    unsigned int input;
    unsigned int extracted;

    friend std::ostream& operator<<(std::ostream& os, const AppLog& log) {
	os << "execution date: " << log.date;
	os << "program version: " << ELASTIC_CLUSTER_SHORT << " " << ELASTIC_CLUSTER_VERSION << "\n";
	os << "program options: " << log.argv << "\n";
	os << "walltime used: " << log.wtime << "\n";
	os << "input edges: " << log.input << "\n";
	os << "extracted clusters: " << log.extracted << "\n";
	return os;
    } // operator<<
}; // struct AppLog


struct edge_t {
    unsigned int id0;
    unsigned int id1;
    unsigned int sim;
}; // struct edge_t

std::ostream& operator<<(std::ostream& os, const edge_t& e) {
    os << "(" << e.id0 << "," << e.id1 << "," << e.sim << ")";
    return os;
} // operator<<

std::istream& operator>>(std::istream& is, edge_t& e) {
    is >> e.id0 >> e.id1 >> e.sim;
    return is;
} // operator>>


void welcome() {
    std::cout << ELASTIC_CLUSTER_FULL << " Version " << ELASTIC_CLUSTER_VERSION << "\n";
    std::cout << ELASTIC_CLUSTER_COPYRIGHT << "\n";
    std::cout << "\n";
} // welcome


std::pair<bool, std::string> run(const AppConfig& opt, AppLog& log, Reporter& report) {
    double t0 = get_time();

    std::ofstream flog((opt.output + ".eclog").c_str());
    if (!flog) return std::make_pair(false, "unable to create " + opt.output + ".eclog");

    std::ofstream fcls((opt.output + ".cls.00000").c_str());
    if (!flog) return std::make_pair(false, "unable to create " + opt.output + ".cls.00000");

    // get files to process
    report << step << "scanning " << opt.input << " for input files..." << std::endl;

    std::vector<fs::path> files;

    if (jaz::files(opt.input, std::back_inserter(files)) == false) {
	return std::make_pair(false, "unable to scan " + opt.input);
    }

    if (files.empty() == true) return std::make_pair(false, "no files to process");

    report << info << "found " << files.size() << " input file(s)" << std::endl;
    report << step << "extracting clusters, be patient..." << std::endl;

    std::vector<unsigned int> uf(opt.nodes);
    unsigned int num_edges = 0;

    jaz::set_make<unsigned int>(&uf[0], uf.size());

    for (unsigned int i = 0; i < files.size(); ++i) {
	if (fs::is_regular_file(files[i]) == false) continue;

	std::ifstream fs(files[i].string().c_str());
	std::istream_iterator<edge_t> ii(fs), end;

	for (; ii != end; ++ii) {
	    if (ii->sim >= opt.level) {
		if (ii->sim > 100) {
		    report.critical << warning << "unexpected edge " << *ii << ", error?" << std::endl;
		}
		if ((ii->id0 < opt.nodes) && (ii->id1 < opt.nodes)) {
		    jaz::set_union(&uf[0], ii->id0, ii->id1);
		    num_edges++;
		} else {
		    report.critical << warning << "incorrect edge " << *ii << ", ignoring" << std::endl;
		}
	    }
	} // for

	if (fs.eof() == false) {
	    report.critical << warning << "incorrect file " << files[i] << ", some edges ignored" << std::endl;
	}

	fs.close();
    } // for i

    std::vector<std::pair<unsigned int, unsigned int> > cluster(opt.nodes);

    for (unsigned int i = 0; i < opt.nodes; ++i) {
	cluster[i].first = jaz::set_find(&uf[0], i);
	cluster[i].second = i;
    }

    std::sort(cluster.begin(), cluster.end());

    unsigned int cnum = 0;
    unsigned int pos = 0;

    for (unsigned int i = 1; i < opt.nodes + 1; ++i) {
	int lim = (i < opt.nodes) ? cluster[i].first : -1;
	if (cluster[pos].first != lim) {
	    for (unsigned int j = pos; j < i - 1; ++j) fcls << cluster[j].second << " ";
	    fcls << cluster[i - 1].second << std::endl;
	    pos = i;
	    cnum++;
	}
    }

    fcls.close();

    log.wtime = get_time() - t0;
    log.input = num_edges;
    log.extracted = cnum;

    report << info << "found " << num_edges << " edges" << std::endl;
    report << info << "extracted " << cnum << " clusters" << std::endl;

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
