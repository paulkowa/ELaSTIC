/***
 *  $Id$
 **
 *  File: disjoint_set.hpp
 *  Created: Oct 04, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2010-2012 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_DISJOINT_SET_HPP
#define JAZ_DISJOINT_SET_HPP

namespace jaz {

  /** Function: set_make
   *  Initializes union-find vector.
   *
   *  Parameters:
   *  uf - array of size n representing union-find vector.
   *  n  - size of union-find vector.
   */
  template <typename Int> inline void set_make(Int* uf, Int n) {
      for (Int i = 0; i < n; ++i) uf[i] = i;
  } // set_make

  /** Function: set_find
   *  Determine which set a particular element is in.
   *
   *  Parameters:
   *  uf - array representing union-find vector.
   *  x  - element to find (must be in the range [0,n), where n is the size of uf).
   */
  template <typename Int> inline Int set_find(Int* uf, Int x) {
      if (uf[x] != uf[uf[x]]) uf[x] = set_find(uf, uf[x]);
      return uf[x];
  } // set_find

  /** Function: set_union
   *  Combine two sets into a single set.
   *
   *  Parameters:
   *  uf - array representing union-find vector.
   *  x  - element from the first set.
   *  y  - element from the second set.
   */
  template <typename Int> inline void set_union(Int* uf, Int x, Int y) {
      Int xx = set_find(uf, x);
      Int yy = set_find(uf, y);
      if (xx != yy ) { if (xx < yy) uf[yy] = xx; else uf[xx] = yy; }
  } // set_union

} // namespace jaz

#endif // JAZ_DISJOINT_SET_HPP
