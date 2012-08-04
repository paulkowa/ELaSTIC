/***
 *  $Id$
 **
 *  File: read_input.hpp
 *  Created: May 23, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the [LICENSE].
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef READ_INPUT_HPP
#define READ_INPUT_HPP

#include <mpi.h>

#include <numeric>

#include <arpa/inet.h>
#include <sys/stat.h>

#include "SequenceCodec.hpp"
#include "Sequence.hpp"
#include "config_log.hpp"


inline unsigned long int file_size(const char* name) {
    struct stat buf;
    int res = stat(name, &buf);
    return (res == 0 ? static_cast<unsigned long int>(buf.st_size) : 0);
} // file_size


inline std::pair<bool, std::string> read_input(const AppConfig& opt, AppLog& log, Reporter& report,
					       MPI_Comm comm, std::vector<Sequence>& seqs) {
    report << step << "reading input sequnces..." << std::endl;

    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // check index size
    unsigned long int fsz = file_size((opt.input + ".eidx").c_str());
    if (fsz == 0) return std::make_pair(false, "unable to open " + opt.input + ".eidx");

    typedef uint16_t index_type;
    std::vector<index_type> index;

    unsigned int n = fsz / sizeof(index_type);
    unsigned int nloc = n / size;

    if (nloc < 2) return std::make_pair(false, "too many processors for this problem");

    // read chunk of the index
    MPI_File fh;
    MPI_Status stat;

    char type[] = "native";
    unsigned int offset = rank * nloc * sizeof(index_type);

    if (rank == size - 1) nloc += n - (nloc * size);
    index.resize(nloc);

    MPI_File_open(comm, const_cast<char*>((opt.input + ".eidx").c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, offset, MPI_BYTE, MPI_BYTE, type, MPI_INFO_NULL);
    MPI_File_read_all(fh, reinterpret_cast<char*>(&index[0]), nloc * sizeof(index_type), MPI_BYTE, &stat);
    MPI_File_close(&fh);

    std::transform(index.begin(), index.end(), index.begin(), ntohs);

    // get offset of the actual data to read
    unsigned int goffset = 0;
    offset = std::accumulate(index.begin(), index.end(), 0);

    MPI_Exscan(&offset, &goffset, 1, MPI_UNSIGNED, MPI_SUM, comm);

    // read actual data
    fsz = file_size((opt.input + ".eseq").c_str());
    if (fsz == 0) return std::make_pair(false, "unable to open " + opt.input + ".eseq");

    std::vector<char> data(offset);

    MPI_File_open(comm, const_cast<char*>((opt.input + ".eseq").c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, goffset, MPI_BYTE, MPI_BYTE, type, MPI_INFO_NULL);
    MPI_File_read_all(fh, &data[0], offset, MPI_BYTE, &stat);
    MPI_File_close(&fh);

    // extract sequences
    SequenceCodec sc;
    seqs.resize(nloc);

    unsigned int pos = index[0];

    seqs[0].id = ntohl(*reinterpret_cast<uint32_t*>(&data[0]));
    seqs[0].s = sc.decode(std::string(data.begin() + 4, data.begin() + pos));

    if ((seqs[0].s.size() < opt.kmer) || (seqs[0].s.find('?') != std::string::npos)) {
	return std::make_pair(false, "invalid input sequence, change kmer?");
    }

    for (unsigned int i = 1; i < nloc; ++i) {
	seqs[i].id = ntohl(*reinterpret_cast<uint32_t*>(&data[pos]));
	seqs[i].s = sc.decode(std::string(data.begin() + pos + 4, data.begin() + pos + index[i]));
	if ((seqs[i].s.size() < opt.kmer) || (seqs[i].s.find('?') != std::string::npos)) {
	    return std::make_pair(false, "invalid input sequence, change kmer?");
	}
	pos += index[i];
    }

    // update log
    log.input = n;

    report << info << "found " << n << " sequences" << std::endl;

    return std::make_pair(true, "");
} // read_input

#endif // READ_INPUT_HPP
