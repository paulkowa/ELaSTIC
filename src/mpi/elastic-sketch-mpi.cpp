/***
 *  $Id$
 **
 *  File: elastic-sketch-mpi.cpp
 *  Created: May 20, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
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

#include "sketch/SequenceDB.hpp"
#include "sketch/config_log.hpp"
#include "sketch/generate_edges.hpp"
#include "sketch/make_shingles.hpp"
#include "sketch/read_input.hpp"
#include "sketch/validate_edges.hpp"

#include <boost/format.hpp>
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
    std::cout << "  --type {nt|aa}        set input sequence type (default nt)\n";
    std::cout << "  --sigma type          use this compressed amino acid alphabet (default A20)\n";
    std::cout << "  --method {0|1}        use this method to validate edges: 0 - alignment, 1 - kmer fraction (default 0)\n";
    std::cout << "  --gaps list           use these parameters for affine gap alignment (default [5,-4,-10,-1])\n";
    std::cout << "  --kmer size           use kmers of this size for sketching (default 15)\n";
    std::cout << "  --level size          use this threshold for edge validation (default 75)\n";
    std::cout << "  --mod size            use this value to perform mod operation in sketching (default 25)\n";
    std::cout << "  --iter size           limit the number of sketching iterations to this size (default 7)\n";
    std::cout << "  --cmax size           use this limit to mark frequent kmers (default 5000)\n";
    std::cout << "  --jmin size           use this limit to extract candidate edges (default 50)\n";
    std::cout << "\n";
} // usage


void run(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm) {
    const int MPI_ABRT_SIG = 13;

    // check MPI comm
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    MPI_Barrier(comm);
    double t0 = MPI_Wtime();

    report << info << "using " << size << " processors..." << std::endl;

    // ALL ERRORS IN THIS FUNCTION ARE CRITICAL AND TERMINATE MPI
    bool res = true;
    std::string err = "";


    // read input sequence
    SequenceList SL;

    boost::tie(res, err) = read_input(opt, log, report, comm, SL);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // initialize RMA
    SequenceRMA rma_seq(comm);
    boost::tie(res, err) = rma_seq.init(SL);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // generate shingles
    jaz::murmur264 hash;
    std::vector<shingle_list_type> shingles;

    boost::tie(res, err) = make_shingles(opt, log, report, SL, hash, shingles);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // generate candidate edges
    std::vector<read_pair> edges;
    boost::tie(res, err) = generate_edges(opt, log, report, comm, SL, shingles, edges);

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // validate candidate edges
    double tv0 = MPI_Wtime() - t0;
    boost::tie(res, err) = validate_edges(opt, log, report, comm, SL, rma_seq, edges);
    double tv1 = MPI_Wtime() - t0;

    if (res == false) {
	report.critical << error << err << std::endl;
	MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // write output
    report << info << "writing output..." << std::endl;

    boost::format fmt("%05d");
    fmt % rank;

    std::ofstream of((opt.output + "." + fmt.str()).c_str());
    std::copy(edges.begin(), edges.end(), std::ostream_iterator<read_pair>(of, "\n"));
    of.close();


    // update log
    double gtv0 = 0.0;
    double gtv1 = 0.0;

    MPI_Reduce(&tv0, &gtv0, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&tv1, &gtv1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    log.through = log.cedges / (gtv1 - gtv0);

    unsigned int l = edges.size();
    MPI_Reduce(&l, &log.vedges, 1, MPI_UNSIGNED, MPI_SUM, 0, comm);

    MPI_Barrier(comm);
    log.wtime = MPI_Wtime() - t0;

    report << info << "generated " << log.vedges << " edges" << std::endl;
    report << info << "edge throughput: " << log.through << std::endl;


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

    report << "time: " << log.wtime << std::endl;
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
