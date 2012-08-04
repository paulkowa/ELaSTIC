/***
 *  $Id$
 **
 *  File: string.hpp
 *  Created: Jun 03, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2012 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_STRING_HPP
#define JAZ_STRING_HPP

#include <string>


namespace jaz {

  /** Function: split
   *  Splits string based on a separator.
   *
   *  Parameters:
   *  pat - Separator character.
   *  s   - String to split.
   *  out - Iterator to store resulting sub-strings.
   */
  template <typename charT, typename traits, typename Alloc, typename Iter>
  void split(charT pat, const std::basic_string<charT, traits, Alloc>& s, Iter out) {
      unsigned int pos = 0;

      for (unsigned int i = 0; i < s.size(); ++i) {
	  if (s[i] == pat) {
	      if (i - pos > 0) {
		  *(out++) = std::basic_string<charT, traits, Alloc>(s, pos, i - pos);
	      }
	      pos = i + 1;
	  }
      } // for

      *(out++) = std::basic_string<charT, traits, Alloc>(s, pos, s.size() - pos);
  } // split


  /** Function: join
   *  Merges a sequence of strings into a single string.
   *  pat   - separator character.
   *  first - Beginning of the sequence to join.
   *  last  - End of the sequence to join.
   *  init  - Prefix to add to the sequence.
   */
  template <typename Iter, typename charT, typename traits, typename Alloc>
  std::basic_string<charT, traits, Alloc>
  join(charT pat, Iter first, Iter last, const std::basic_string<charT, traits, Alloc>& init) {
      std::basic_string<charT, traits, Alloc> s(init);
      if (s.empty() == true) s = *(first++);
      for (; first != last; ++first) s += pat + std::basic_string<charT, traits, Alloc>(*first);
      return s;
  } // join

  /** Function: join
   */
  template <typename Iter> inline std::string join(char pat, Iter first, Iter last) {
      return join(pat, first, last, std::string(""));
  } // join


  /** Function: approx_match
   *  Performs approximate string matching.
   *
   *  Parameters:
   *  T  - Text string.
   *  P  - Pattern string.
   *  mm - Allowed mismatches.
   *
   *  Returns:
   *  Pair in which the first element stores position of the first match or
   *  std::string::npos if no match is found, and the second element is
   *  true if more than one match exists.
   */
  template <typename charT, typename traits, typename Alloc>
  std::pair<long int, bool> approx_match(const std::basic_string<charT, traits, Alloc>& T,
					 const std::basic_string<charT, traits, Alloc>& P,
					 unsigned int mm) {
      unsigned int n = T.size();
      unsigned int m = P.size();

      long int pos = std::string::npos;
      bool mult = false;

      long int l = n - m + 1;
      unsigned int cmm = mm + 1;

      for (long int i = 0; i < l; ++i) {
	  unsigned int t = 0;
	  unsigned int j = 0;

	  for (; j < m; ++j) {
	      if (T[i + j] != P[j]) t++;
	      if (t == cmm) break;
	  }

	  if (j == m) {
	      if (t < cmm - 1) {
		  pos = i;
		  mult = false;
		  cmm = t + 1;
	      } else if (t == cmm - 1) {
		  if (pos == std::string::npos) pos = i;
		  else mult = true;
	      }
	  }
      } // for i

      return std::make_pair(pos, mult);
  } // approx_match

} // namespace jaz

#endif // JAZ_STRING_HPP
