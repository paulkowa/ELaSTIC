/***
 *  $Id$
 **
 *  File: files.hpp
 *  Created: Mar 26, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_FILES_HPP
#define JAZ_FILES_HPP

#include <algorithm>
#include <string>

#include <boost/filesystem.hpp>


namespace fs = boost::filesystem;

namespace jaz {

  template <typename OutputIter>
  inline bool files(const std::string& dir, OutputIter out) {
      fs::path p(dir);
      if (fs::exists(p) == false) return false;

      if (fs::is_directory(p) == true) {
	  std::copy(fs::directory_iterator(p),fs::directory_iterator(), out);
      } else (*out) = p;

      return true;
  } // files

} // namespace jaz

#endif // JAZ_FILES_HPP
