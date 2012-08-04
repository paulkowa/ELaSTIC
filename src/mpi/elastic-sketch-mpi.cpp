/***
 *  $Id$
 **
 *  File: elastic-sketch-mpi.cpp
 *  Created: May 20, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#include <mpi.h>

#include <iostream>
#include <map>
#include <string>

#include "config.hpp"
#include "iomanip.hpp"

#include "sketch/Sequence.hpp"
#include "sketch/config_log.hpp"
#include "sketch/generate_edges.hpp"
#include "sketch/make_shingles.hpp"
#include "sketch/read_input.hpp"
#include "sketch/validate_edges.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <jaz/parameters.hpp>
#include <jaz/hash.hpp>


void welcome() {
    std::cout << ELASTIC_SKETCH_FULL << " Version " << ELASTIC_SKETCH_VERSION << "\n";
    std::cout << ELASTIC_SKETCH_COPYRIGHT << "\n";
    std::cout << std::endl;
} // welcome

void usage() {
    std::cout << "Usage: elastic-sketch-mpi --input name --output name [options...]\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  --input name          read input from files with this prefix\n";
    std::cout << "  --output name         write output to this file\n";
    std::cout << "  --config name         read configuration from this file\n";
    std::cout << "  --kmer size           use kmers of this size for sketching (default 15)\n";
    std::cout << "  --mod size            use this value to perform mod operation in sketching (default 25)\n";
    std::cout << "  --iter size           limit the number of sketching iterations to this size (default -1)\n";
    std::cout << "  --cmax size           use this limit to mark frequent kmers (default 5000)\n";
    std::cout << "  --jmin size           use this limit to filter similar sequences (default 50)\n";
    std::cout << "\n";
} // usage


void run(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm) {
    const int MPI_ABRT_SIG = 13;

    // check MPI comm
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // ALL ERRORS IN THIS FUNCTION ARE CRITICAL AND TERMINATE MPI
    bool res = true;
    std::string err = "";


    // read input sequence
    std::vector<Sequence> seqs;
    boost::tie(res, err) = read_input(opt, log, report, comm, seqs);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // generate shingles
    jaz::murmur264 hash;
    std::vector<shingle_list_type> shingles;

    boost::tie(res, err) = make_shingles(opt, log, report, seqs, hash, shingles);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // generate candidate edges
    std::vector<read_pair> edges;
    boost::tie(res, err) = generate_edges(opt, log, report, seqs, shingles, comm, edges);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // validate candidate edges
    boost::tie(res, err) = validate_edges(opt, log, report, seqs, edges, comm);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // done :-)

    std::ofstream of((opt.output + "." + boost::lexical_cast<std::string>(rank)).c_str());
    std::copy(edges.begin(), edges.end(), std::ostream_iterator<read_pair>(of, "\n"));
    of.close();

    // write log
    if (rank == 0) {
	std::ofstream flog((opt.output + ".eslog").c_str());

	if (!flog) {
	    report.critical << error << "unable to create " << opt.output + ".eslog" << std::endl;
	    MPI_Abort(comm, MPI_ABRT_SIG);
	}

	flog << log;
	flog.close();
    }

    MPI_Barrier(comm);

    report << "done!" << std::endl;
} // run


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) welcome();

    // get parameters
    std::map<std::string, std::string> conf;

    if (argc == 1) {
	if (rank == 0) usage();
	return MPI_Finalize();
    }

    bool res = false;
    int pos = 0;

    boost::tie(res, pos) = jaz::parse_argv(argc, argv, conf);

    if (res == false) {
	if (pos == -1) {
	    if (rank == 0) {
		usage();
		std::cout << error << "incorrect command line arguments\n";
	    }
	    return MPI_Finalize();
	} else {
	    if (rank == 0) {
		std::cout << error << "incorrect command line argument " << argv[pos] << "\n";
	    }
	    return MPI_Finalize();
	}
    } // if res

    // set reporter
    Reporter report = (rank == 0) ? Reporter(std::cout, std::cout) : Reporter(nout, std::cout);

    // create config and log
    AppConfig opt;
    std::string err = "";

    boost::tie(res, err) = opt.set(conf);

    if (res == false) {
	report << error << err << "\n";
	return MPI_Finalize();
    }

    AppLog log;
    log.argv = jaz::join(' ', argv + 1, argv + argc);
    log.cpus = size;

    // and here we go
    run(opt, log, report, MPI_COMM_WORLD);

    return MPI_Finalize();
} // main
