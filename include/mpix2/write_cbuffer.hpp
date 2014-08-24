/***
 *  $Id$
 **
 *  File: write_cbuffer.hpp
 *  Created: Feb 12, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of mpix2.
 */

#ifndef MPIX2_WRITE_CBUFFER_HPP
#define MPIX2_WRITE_CBUFFER_HPP

#include <string>
#include <mpi.h>


namespace mpix {

  inline bool write_cbuffer(const std::string& name, const char* buf, int bsz,
			    MPI_Comm comm) {
      MPI_File fh;
      MPI_Status stat;

      int rank;
      MPI_Comm_rank(comm, &rank);

      if (rank == 0) MPI_File_delete(const_cast<char*>(name.c_str()), MPI_INFO_NULL);
      MPI_Barrier(comm);

      int err = MPI_File_open(comm, const_cast<char*>(name.c_str()),
			      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

      if (err != MPI_SUCCESS) return false;

      unsigned long long int size = bsz;
      unsigned long long int offset = 0;

      MPI_Scan(&size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
      offset -= size;

      char type[] = "native";

      MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, type, MPI_INFO_NULL);
      MPI_File_write_all(fh, const_cast<char*>(buf), size, MPI_CHAR, &stat);
      MPI_File_close(&fh);

      return true;
  } // write_cbuffer

} // namespace mpix

#endif // MPIX2_WRITE_CBUFFER_HPP
