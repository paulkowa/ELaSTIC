/***
 *  $Id$
 **
 *  File: parameters.hpp
 *  Created: Dec 12, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2007-2012 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of jaz.
 */

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <algorithm>
#include <cctype>
#include <fstream>
#include <string>
#include <vector>
#include "string.hpp"


namespace jaz {

  /** Function: read_config
   */
  template <typename Container>
  std::pair<bool, int> read_config(const std::string& name, Container& cfg) {
      std::ifstream f(name.c_str());
      if (!f) return std::make_pair(false, -1);

      int line = 0;
      std::vector<std::string> tokens;

      do {
          std::string s;
          std::getline(f, s);
          ++line;

          if ((f.good() == false) && (f.eof() == false)) return std::make_pair(false, line);

          if (s.empty() == false) {
              // needed for windows files
              unsigned int len = s.size() - 1;
              if (s[len] == '\r') s.resize(len);

              if (s.empty() == true) continue;

              if ((s[0] != '#') && (s[0] != ';')) {
                  tokens.clear();
                  split('=', s, std::back_inserter(tokens));
                  if (tokens.size() != 2) return std::make_pair(false, line);
                  else {
                      tokens[0].erase(std::remove_if(tokens[0].begin(), tokens[0].end(), isspace), tokens[0].end());
                      tokens[1].erase(std::remove_if(tokens[1].begin(), tokens[1].end(), isspace), tokens[1].end());
                      cfg[tokens[0]] = tokens[1];
                  }
              }
          } // if s.empty()
      }
      while (f.eof() == false);

      f.close();
      return std::make_pair(true, line);
  } // read_config

  /** Function: parse_argv
   */
  template <typename Container>
  std::pair<bool, int> parse_argv(int argc, char* argv[], Container& cfg) {
      if (argc % 2 != 1) return std::make_pair(false, -1);

      int pos = 1;

      for (; pos < argc; pos += 2) {
          std::string opt = argv[pos];
          std::string arg = argv[pos + 1];
          if ((opt.size() < 3) || (opt[0] != '-') || (opt[1] != '-')) return std::make_pair(false, pos);
          cfg[std::string(opt.begin() + 2, opt.end())] = arg;
      } // for pos

      return std::make_pair(true, pos);
  } // parse_argv

  /** Function: check_option
   */
  template <typename Key, typename Val, typename Container>
  inline bool check_option(const Container& cfg, const Key& key, Val& val) {
      typename Container::const_iterator iter(cfg.find(key));
      if (iter == cfg.end()) return false;
      val = iter->second;
      return true;
  } // check_option

} // namespace jaz

#endif // PARAMETERS_HPP
