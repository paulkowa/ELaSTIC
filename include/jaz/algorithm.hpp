/***
 *  $Id$
 **
 *  File: algorithm.hpp
 *  Created: Sep 09, 2011
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_ALGORITHM_HPP
#define JAZ_ALGORITHM_HPP

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <functional>
#include <iterator>
#include <numeric>


namespace jaz {

  /** Function: copy_n
   */
  template <typename InputIter, typename Size, typename OutputIter>
  OutputIter copy_n(InputIter first, InputIter last, Size n, OutputIter out) {
      if (n > 0) {
	  Size count = 0;
	  for (; (first != last) && (count < n); ++first, ++count) {
	      *out++ = *first;
	  } // for
      } // if
      return out;
  } // copy_n


  /** Function: compact
   */
  template <typename Iter, typename Oper>
  Iter compact(Iter first, Iter last, Oper op) {
      if (first == last) return last;

      Iter res = first;
      Iter pos = first;

      first++;

      for (; first != last; ++first) {
	  if (*pos < *first) {
	      *res = *pos;
	      pos++;
	      *res = std::accumulate(pos, first, *res, op);
	      pos = first;
	      res++;
	  }
      }

      *res = *pos;
      pos++;
      *res = std::accumulate(pos, last, *res, op);
      res++;

      return res;
  } // compact


  /** Function: mode
   */
  template <typename Iter, typename Pred>
  std::pair<Iter, Iter> mode(Iter first, Iter last, Pred pred) {
      if (first == last) return std::make_pair(first, first);

      std::size_t count = 0;
      Iter res_first;
      Iter res_last;

      std::size_t cur_count = 1;
      Iter cur_first = first;

      first++;

      for (; first != last; ++first) {
	  if (pred(*cur_first, *first)) {
	      if (count < cur_count) {
		  res_first = cur_first;
		  res_last = first;
		  count = cur_count;
	      }
	      cur_first = first;
	      cur_count = 1;
	  } else cur_count++;
      }

      if (count < cur_count) {
	  res_first = cur_first;
	  res_last = first;
      }

      return std::make_pair(res_first, res_last);
  } // mode

  template <typename Iter>
  std::pair<Iter, Iter> mode(Iter first, Iter last) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return mode(first, last, std::less<value_type>());
  } // mode


  /** Function: intersection_size
   */
  template <typename Iter1, typename Iter2, typename Pred>
  std::size_t intersection_size(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2, Pred pred) {
      std::size_t S = 0;

      while ((first1 != last1) && (first2 != last2)) {
	  if (pred(*first1, *first2)) ++first1;
	  else if (pred(*first2, *first1)) ++first2;
	  else {
	      first1++;
	      first2++;
	      S++;
	  }
      } // while

      return S;
  } // intersection_size

  /** Function: intersection_size
   */
  template <typename Iter1, typename Iter2>
  inline std::size_t intersection_size(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) {
      typedef typename std::iterator_traits<Iter1>::value_type value_type;
      return intersection_size(first1, last1, first2, last2, std::less<value_type>());
  } // intersection_size


  /** Function: count_unique
   */
  template <typename Iter, typename Pred>
  std::size_t count_unique(Iter first, Iter last, Pred pred) {
      if (first == last) return 0;
      std::size_t S = 1;
      Iter prev = first++;
      for (; first != last; ++first, ++prev) if (pred(*prev, *first) == false) ++S;
      return S;
  } // count_unique

  /** Function: count_unique
   */
  template <typename Iter>
  inline std::size_t count_unique(Iter first, Iter last) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return count_unique(first, last, std::equal_to<value_type>());
  } // count_unique


  /** Function: find_all
   */
  template <typename InputIter, typename OutputIter, typename Pred>
  OutputIter find_all(InputIter first, InputIter last, OutputIter out, Pred pred) {
      for (; first != last; ++first) {
	  if (pred(*first) == true) {
	      *out = first;
	      ++out;
	  }
      } // for
      return out;
  } // find_all


  /** Function: random_sample_n
   */
  template <typename InputIter, typename OutputIter, typename Size, typename Random>
  OutputIter random_sample_n(InputIter first, InputIter last, OutputIter out, Size n, Random& rand) {
      Size m = std::distance(first, last);
      Size k = std::min(n, m);

      while (k > 0) {
  	  if (rand() < k) {
  	      *out = *first;
  	      ++out;
  	      --k;
  	  }
  	  --m;
  	  ++first;
      } // while

      return out;
  } // random_sample_n


  /** Function: max_vectors
   */
  template <typename Iter, typename Pred>
  Iter max_vectors(Iter first, Iter last, Pred pred) {
      if (first == last) return first;

      Iter beg = first;
      Iter end = last;

      do {
	  Iter cur = beg;
	  Iter i = beg;

	  i++;

	  // find max element
	  for (; i != end; ++i) {
	      if (pred(*i, *beg) == true) {
		  --end; std::swap(*i, *end);
		  --i;
	      } else if (pred(*beg, *i) == true) {
		  std::swap(*i, *beg);
		  cur = i;
		  cur++;
	      }
	  } // for i

	  i = beg;
	  i++;

	  // clean
	  for (; (i != cur) && (i != end); ++i) {
	      if (pred(*i, *beg) == true) {
		  --end; std::swap(*i, *end);
		  --i;
	      }
	  } // for i

	  beg++;
      } while (beg != end);

      return end;
  } // max_vectors

} // namespace jaz

#endif // JAZ_ALGORITHM_HPP
