/***
 *  $Id$
 **
 *  File: elastic-sketch.cpp
 *  Created: May 20, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#include <mpi.h>

#include <iostream>
#include <map>
#include <new>
#include <string>
#include <sstream>

#include "config.hpp"
#include "iomanip.hpp"

#include "sketch/SequenceDB.hpp"
#include "sketch/config_log.hpp"
#include "sketch/generate_edges.hpp"
#include "sketch/make_shingles.hpp"
#include "sketch/read_input.hpp"
#include "sketch/validate_edges.hpp"
#include "sketch/write_output.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
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


void run(const AppConfig& opt, AppLog& log, Reporter& report, MPI_Comm comm) {
    const int MPI_ABRT_SIG = 13;

    // check MPI comm
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    MPI_Barrier(comm);
    double t0 = MPI_Wtime();

    report << info << "using " << size << " processors..." << std::endl;
    report << info << "input configuration...\n" << opt;

    // ALL ERRORS IN THIS FUNCTION ARE CRITICAL AND TERMINATE MPI
    bool res = true;
    std::string err = "";


    // set random seed
    std::srand((rank + 1) * std::time(0));


    // read input sequence
    SequenceList SL;

    boost::tie(res, err) = read_input(opt, log, report, comm, SL);

    if (res == false) {
        report.critical << error << err << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }

    // generate profiling report
    if (rank == opt.dbg) {
        using namespace boost::accumulators;
        accumulator_set<unsigned int, stats<tag::min, tag::max, tag::mean, tag::sum> > acc;

        for (unsigned int i = 0; i < SL.seqs.size(); ++i) acc(SL.seqs[i].s.size());
        unsigned int mem = SL.seqs.size() * sizeof(Sequence) + sum(acc);

        report.stream << debug << "SL.seqs: " << SL.seqs.size()
                      << ", min/avg/max: " << min(acc) << "/" << mean(acc) << "/" << max(acc)
                      << " memory: " << to_size(mem) << std::endl;
    } // if dbg


    // generate shingles
    jaz::murmur264 hash;
    std::vector<shingle_list_type> shingles;

    boost::tie(res, err) = make_shingles(opt, log, report, SL, hash, shingles);

    if (res == false) {
        report.critical << error << err << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }

    // generate profiling report
    if (rank == opt.dbg) {
        unsigned int S = 0;

        for (unsigned int i = 0; i < shingles.size(); ++i) S += shingles[i].size();
        unsigned int mem = sizeof(shingles) + S * sizeof(shingle_list_type::value_type);

        report.stream << debug << "shingles, memory: " << to_size(mem) << std::endl;
    } // if dbg


    // compress sequences
    if ((opt.is_dna == false) && (opt.compress == true) && (opt.sigma != "A20")) {
        CompressedAlphabet ca(opt.sigma);
        unsigned int n = SL.seqs.size();
        for (unsigned int i = 0; i < n; ++i) {
            SL.seqs[i].s = ca(SL.seqs[i].s);
        }
    }


    // generate candidate edges
    // at the end processors are most likely out of sync
    std::vector<read_pair> edges;

    boost::tie(res, err) = generate_edges(opt, log, report, comm, SL, shingles, edges);

    if (res == false) {
        report.critical << error << err << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }

    unsigned long long int etot = edges.size();
    unsigned long long int num_cedges = 0;

    MPI_Reduce(&etot, &num_cedges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);

    report << info << "extracted " << num_cedges << " candidate edges" << std::endl;

    if (rank == opt.dbg) {
        report.stream << debug << "collective edge memory: " << to_size(num_cedges * sizeof(read_pair)) << std::endl;
    }

#ifdef WITH_MPE
    // profiling of validation is disabled by default
    return;
#endif

    // initialize RMA
    SequenceRMA rma_seq(comm);
    boost::tie(res, err) = rma_seq.init(SL);

    if (res == false) {
        report.critical << error << err << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // validate candidate edges
    double tv0 = MPI_Wtime() - t0;
    if (opt.validate == true) {
        if (opt.wsq == false) boost::tie(res, err) = validate_edges(opt, log, report, comm, SL, rma_seq, edges);
        else boost::tie(res, err) = validate_edges_ws(opt, log, report, comm, SL, rma_seq, edges);
    } else {
        report << info << "skipping validation..." << std::endl;
    }
    double tv1 = MPI_Wtime() - t0;

    if (res == false) {
        report.critical << error << err << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // this barrier is needed for BG
    report << info << "syncing..." << std::endl;
    MPI_Barrier(comm);


    // write output
    double ti0 = MPI_Wtime();
    unsigned long int fs = 0;
    boost::tie(res, err) = write_output(opt, log, report, comm, edges, fs);
    double ti1 = MPI_Wtime();

    if (res == false) {
        report.critical << error << err << std::endl;
        MPI_Abort(comm, MPI_ABRT_SIG);
    }


    // final sync
    MPI_Barrier(comm);
    log.wtime = MPI_Wtime() - t0;


    // update log
    double gtv0 = 0.0;
    double gtv1 = 0.0;

    MPI_Reduce(&tv0, &gtv0, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&tv1, &gtv1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    double gst = 0.0;
    MPI_Reduce(&tv0, &gst, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    log.cedges = num_cedges;

    unsigned long long int l = edges.size();
    MPI_Reduce(&l, &log.vedges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);

    double fa2a = (2.0 * log.cedges) / ((1.0 * SL.N) * (SL.N - 1));

    report << info << "extracted " << log.cedges << " candidate edges" << std::endl;
    report << info << "generated " << log.vedges << " valid edges" << std::endl;
    report << info << "sketching phase: " << gst << std::endl;

    if (opt.validate == true) {
        report << info << "validation phase: " << (gtv1 - gtv0) << std::endl;
        report << info << "edge throughput: " << static_cast<unsigned int>(num_cedges / (gtv1 - gtv0)) << std::endl;
        report << info << "fraction processed: " << fa2a << std::endl;
        report << info << "fraction accepted: " << static_cast<double>(log.vedges) / log.cedges << std::endl;
    }

    // write log
    if (rank == 0) {
        std::ofstream flog((opt.output + ".eslog").c_str());

        if (!flog) {
            report.critical << error << "unable to create " << opt.output + ".eslog" << std::endl;
            MPI_Abort(comm, MPI_ABRT_SIG);
        }

        flog << log;
        flog << "config:" << std::endl;
        flog << opt;
        flog.close();
    }

    report << "time: " << log.wtime << " [s]" << std::endl;
    report << "i/o time: " << (ti1 - ti0) << " [s]" << std::endl;
    report << "i/o band: " << fs / ((ti1 - ti0) * 1024 * 1024) << " [MB/s]" << std::endl;
    report << "su: " << (size * log.wtime) / 3600 << std::endl;
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
        if (rank == 0) AppConfig::usage();
        return MPI_Finalize();
    }

    bool res = false;
    int pos = 0;

    boost::tie(res, pos) = jaz::parse_argv(argc, argv, conf);

    if (res == false) {
        if (pos == -1) {
            if (rank == 0) {
                AppConfig::usage();
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
    Reporter report = (rank == 0) ? Reporter(std::cout, std::cout) : Reporter(std::cout, nout, std::cout);

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
