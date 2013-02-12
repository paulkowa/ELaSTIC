/***
 *  $Id$
 **
 *  File: write_cbuffer.hpp
 *  Created: Feb 12, 2013
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of mpix.
 */

#ifndef MPIX_WRITE_CBUFFER_HPP
#define MPIX_WRITE_CBUFFER_HPP

#include <string>
#include <mpi.h>


namespace mpix {

  inline bool write_cbuffer(const std::string& name, const char* buf, unsigned int bsz, MPI_Comm comm) {
      MPI_File fh;
      MPI_Status stat;

      int err = MPI_File_open(comm, const_cast<char*>(name.c_str()),
			      MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

      if (err != MPI_SUCCESS) return false;

      unsigned int size = bsz;
      unsigned int offset = 0;

      MPI_Scan(&size, &offset, 1, MPI_UNSIGNED, MPI_SUM, comm);
      offset -= size;

      char type[] = "native";

      MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, type, MPI_INFO_NULL);
      MPI_File_write_all(fh, const_cast<char*>(buf), size, MPI_CHAR, &stat);
      MPI_File_close(&fh);

      return true;
  } // write_cbuffer

} // namespace mpix

#endif // MPIX_WRITE_CBUFFER_HPP
