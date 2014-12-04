/***
 *  $Id$
 **
 *  File: algorithm.hpp
 *  Created: Sep 09, 2011
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2004-2013 Jaroslaw Zola
 *  Distributed under the Boost Software License, Version 1.0.
 *  See accompanying file LICENSE_BOOST.txt.
 *
 *  This file is part of jaz0x.
 */

#ifndef JAZ_ALGORITHM_HPP
#define JAZ_ALGORITHM_HPP

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
// #include <random>


namespace jaz {

  /** Function: copy_n
   *  This function is different from std::copy_n.
   *  It copies at most n elements!
   */
  template <typename InputIter, typename Size, typename OutputIter>
  OutputIter copy_n(InputIter first, InputIter last, Size n, OutputIter out) {
      if (n > 0) {
          for (; (first != last) && (n > 0); ++first, --n) *out++ = *first;
      }
      return out;
  } // copy_n


  /** Function: count_unique
   */
  template <typename Iter, typename Comp>
  std::size_t count_unique(Iter first, Iter last, Comp comp) {
      if (first == last) return 0;
      std::size_t S = 1;
      Iter prev = first++;
      for (; first != last; ++first, ++prev) if (comp(*prev, *first) == false) ++S;
      return S;
  } // count_unique

  /** Function: count_unique
   */
  template <typename Iter>
  inline std::size_t count_unique(Iter first, Iter last) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return count_unique(first, last, std::equal_to<value_type>());
  } // count_unique


  /** Function: range
   */
  template <typename Iter, typename Comp>
  Iter range(Iter first, Iter last, Comp comp) {
      if (first == last) return last;
      Iter iter = first;
      for (; (iter != last) && comp(*first, *iter); ++iter);
      return iter;
  } // range

  template <typename Iter>
  inline Iter range(Iter first, Iter last) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return range(first, last, std::equal_to<value_type>());
  } // range


  /** Function: compact
   */
  template <typename Iter, typename Oper, typename Comp>
  Iter compact(Iter first, Iter last, Oper op, Comp comp) {
      if (first == last) return last;

      Iter res = first;

      while (first != last) {
          Iter iter = range(first, last, comp);
          *res = *first;
          ++first;
          *res = std::accumulate(first, iter, *res, op);
          first = iter;
          ++res;
      }

      return res;
  } // compact

  template <typename Iter, typename Oper>
  inline Iter compact(Iter first, Iter last, Oper op) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return compact(first, last, op, std::equal_to<value_type>());
  } // compact


  /** Function: mode
   *  Finds the first range containing mode.
   */
  template <typename Iter, typename Comp>
  std::pair<Iter, Iter> mode(Iter first, Iter last, Comp comp) {
      if (first == last) return std::make_pair(last, last);

      std::pair<Iter, Iter> res = std::make_pair(last, last);
      std::size_t mcount = 1;

      Iter iter = first;

      while (first != last) {
          std::size_t count = 0;
          for (; (iter != last) && comp(*first, *iter); ++iter, ++count);
          if (mcount < count) {
              res.first = first;
              res.second = iter;
              mcount = count;
          }
          first = iter;
      }

      return res;
  } // mode

  template <typename Iter>
  inline std::pair<Iter, Iter> mode(Iter first, Iter last) {
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      return mode(first, last, std::equal_to<value_type>());
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
              ++first1;
              ++first2;
              ++S;
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
  /*
  template <typename InputIter, typename OutputIter, typename Size, typename Random>
  OutputIter random_sample_n(InputIter first, InputIter last, OutputIter out, Size n, Random&& rand) {
      Size m = std::distance(first, last);
      Size k = std::min(n, m);

      while (k > 0) {
          if (rand() % m < k) {
              *out = *first;
              ++out;
              --k;
          }
          --m;
          ++first;
      } // while

      return out;
  } // random_sample_n
  */

  /** Function: random_sample_n
   */
  /*
  template <typename InputIter, typename OutputIter, typename Size>
  inline OutputIter random_sample_n(InputIter first, InputIter last, OutputIter out, Size n) {
      std::random_device rd;
      std::mt19937 rand{rd()};
      return random_sample_n(first, last, out, n, rand);
  } // random_sample_n
  */


  /** Class: reservoir_sampler
   *  Functor to sample k elements from a stream of unknown length.
   */
  /*
  template <typename T, typename Sequence = std::vector<T>, typename Random = std::mt19937>
  class reservoir_sampler {
  public:
      reservoir_sampler(int k, Sequence& sample) : k_(k), sample_(sample), count_(0), rand_(lrand_) {
          first_ = sample_.size();
          std::random_device rd;
          lrand_ = Random(rd());
      } // reservoir_sampler

      reservoir_sampler(int k, Sequence& sample, Random& rand)
          : reservoir_sampler(k, sample) { rand_ = rand; }

      void operator()(const T& t) {
          if (count_ < k_) sample_.push_back(t);
          else {
              udist_ = std::uniform_int_distribution<long int>{0, count_ - 1};
              int pos = udist_(rand_);
              if (pos < k_) sample_[first_ + pos] = t;
          }
          ++count_;
      } // operator()

      template <typename Iter>
      void operator()(Iter first, Iter last) {
          for (auto i = first; i != last; ++i) this->operator()(*i);
      } // operator()

  private:
      int k_;
      int first_;
      Sequence& sample_;

      long int count_;

      Random& rand_;
      Random lrand_;
      std::uniform_int_distribution<long int> udist_;

  }; // class reservoir_sampler

  template <typename Sequence>
  reservoir_sampler<typename Sequence::value_type, Sequence> make_reservoir_sampler(int k, Sequence& sample) {
      return reservoir_sampler<typename Sequence::value_type, Sequence>(k, sample);
  } // make_reservoir_sampler

  template <typename Sequence, typename Random>
  reservoir_sampler<typename Sequence::value_type, Sequence, Random> make_reservoir_sampler(int k, Sequence& sample, Random& rand) {
      return reservoir_sampler<typename Sequence::value_type, Sequence, Random>(k, sample, rand);
  } // make_reservoir_sampler
  */


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

          ++i;

          // find max element
          for (; i != end; ++i) {
              if (pred(*i, *beg) == true) {
                  --end; std::swap(*i, *end);
                  --i;
              } else if (pred(*beg, *i) == true) {
                  std::swap(*i, *beg);
                  cur = i;
                  ++cur;
              }
          } // for i

          i = beg;
          ++i;

          // clean
          for (; (i != cur) && (i != end); ++i) {
              if (pred(*i, *beg) == true) {
                  --end; std::swap(*i, *end);
                  --i;
              }
          } // for i

          ++beg;
      } while (beg != end);

      return end;
  } // max_vectors

} // namespace jaz

#endif // JAZ_ALGORITHM_HPP
