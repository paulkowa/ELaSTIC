/***
 *  $Id$
 **
 *  File: sequence_compare.hpp
 *  Created: May 03, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *
 *  Boost Software License - Version 1.0 - August 17th, 2003
 *
 *  Permission is hereby granted, free of charge, to any person or organization
 *  obtaining a copy of the software and accompanying documentation covered by
 *  this license (the "Software") to use, reproduce, display, distribute,
 *  execute, and transmit the Software, and to prepare derivative works of the
 *  Software, and to permit third-parties to whom the Software is furnished to
 *  do so, all subject to the following:
 *
 *  The copyright notices in the Software and this entire statement, including
 *  the above license grant, this restriction and the following disclaimer,
 *  must be included in all copies of the Software, in whole or in part, and
 *  all derivative works of the Software, unless such copies or derivative
 *  works are solely in the form of machine-executable object code generated by
 *  a source language processor.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 *  SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 *  FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 *  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 */

#ifndef SEQUENCE_COMPARE_HPP
#define SEQUENCE_COMPARE_HPP

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>


namespace bio {

  namespace detail {

    // this code comes from jaz
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


    template <typename Iter1, typename Iter2>
    int count_distance(Iter1 first1, Iter1 last1, Iter2 first2, Iter2 last2) {
	int S = 0;

	while ((first1 != last1) && (first2 != last2)) {
	    if (first1->first < first2->first) {
		S += (first1->second * first1->second);
		++first1;
	    }
	    else if (first2->first < first1->first) {
		S += (first2->second * first2->second);
		++first2;
	    }
	    else {
		int d = (first1->second - first2->second);
		S += d * d;
		first1++;
		first2++;
	    }
	} // while

	return S;
    } // count_distance


    template <typename Sequence> void general_kmer_index(const std::string& s, unsigned int k, Sequence& S) {
	unsigned int l = s.size();
	unsigned int end = l - k + 1;
	S.resize(end);
	for (unsigned int i = 0; i < end; ++i) {
	    S[i] = std::string(s.begin() + i, s.begin() + i + k);
	}
	std::sort(S.begin(), S.end());
    } // general_kmer_index


    template <typename Map> void general_kmer_count(const std::string& s, unsigned int k, Map& S) {
	S.clear();
	unsigned int l = s.size();
	unsigned int end = l - k + 1;
	for (unsigned int i = 0; i < end; ++i) {
	    S[std::string(s.begin() + i, s.begin() + i + k)]++;
	}
    } // general_kmer_count


    class dna_digit {
    public:
	dna_digit() {
	    std::memset(digit_, 0, 256);
	    digit_['c'] = digit_['C'] = 1;
	    digit_['g'] = digit_['G'] = 2;
	    digit_['t'] = digit_['T'] = 3;
	    digit_['u'] = digit_['U'] = 3;
	} // dna_digit

    protected:
	char digit_[256];

    }; // dna_digit


    class dna_kmer_index : public dna_digit {
    public:
	dna_kmer_index() : dna_digit() { }

	template <typename Sequence>
	void operator()(const std::string& s, unsigned int k, Sequence& S) {
	    unsigned int l = s.size();
	    unsigned int end = l - k + 1;

	    S.resize(end);

	    // first kmer
	    unsigned long long int v = digit_[s[k - 1]];
	    for (unsigned int i = 0; i < k - 1; ++i) {
		v += digit_[s[i]] * (1ULL << ((k - i - 1) << 1));
	    }

	    S[0] = v;

	    // and then all other
	    unsigned long long int b = 1ULL << ((k - 1) << 1);

	    for (unsigned int i = 1; i < end; ++i) {
		v = (v - b * digit_[s[i - 1]]) * 4 + digit_[s[i + k - 1]];
		S[i] = v;
	    }

	    std::sort(S.begin(), S.end());
	} // operator()

    }; // class dna_kmer_index


    class dna_kmer_count : public dna_digit {
    public:
	dna_kmer_count() : dna_digit() { }

	template <typename Map>
	void operator()(const std::string& s, unsigned int k, Map& S) {
	    unsigned int l = s.size();
	    unsigned int end = l - k + 1;

	    // first kmer
	    unsigned long long int v = digit_[s[k - 1]];

	    for (unsigned int i = 0; i < k - 1; ++i) {
		v += digit_[s[i]] * (1ULL << ((k - i - 1) << 1));
	    }

	    S[v] = 1;

	    // and then all other
	    unsigned long long int b = 1ULL << ((k - 1) << 1);

	    for (unsigned int i = 1; i < end; ++i) {
		v = (v - b * digit_[s[i - 1]]) * 4 + digit_[s[i + k - 1]];
		S[v]++;
	    }
	} // operator()

    }; // class dna_kmer_count

  } // namespace detail



  template <typename Derived> struct sequence_compare {
      boost::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
	  return static_cast<Derived*>(this)->operator()(s0, s1);
      }
  }; // struct sequence_compare



  /** Class: scoring_matrix
   *
   *  Functor encapsulating a scoring_matrix functionality.
   */
  class scoring_matrix {
  public:
      scoring_matrix() : sz_(0), matrix_(0) { std::memset(sigma_, 0, 256); }

      /** Constructor: scoring_matrix
       *
       *  Parameter:
       *  sigma -  Map of the alphabet used by the matrix.
       *  matrix - Raw-wise substitution matrix.
       */
      scoring_matrix(unsigned char sigma[256], const std::vector<char>& matrix)
	  : matrix_(matrix), sz_(std::sqrt(matrix.size())) { std::memcpy(sigma_, sigma, 256); }

      /** Function: operator()
       *
       *  Returns:
       *  Substitution score between a and b.
       */
      int operator()(char a, char b) const { return matrix_[sigma_[a] * sz_ + sigma_[b]]; }


  private:
      unsigned int sz_;
      unsigned char sigma_[256];
      std::vector<char> matrix_;

  }; // scoring_matrix


  scoring_matrix make_dummy_sm(int m, int s) {
      unsigned char sigma[256];
      for (unsigned int i = 0; i < 256; ++i) sigma[i] = i;
      std::vector<char> matrix(256 * 256, s);
      for (unsigned int i = 0; i < 256; ++i) matrix[(i << 8) + i] = m;
      return scoring_matrix(sigma, matrix);
  } // make_dummy_sm

  scoring_matrix make_dna_sm(int m, int s) {
      unsigned char sigma[256];
      for (unsigned int i = 0; i < 256; ++i) sigma[i] = 4;

      sigma['a'] = sigma['A'] = 0;
      sigma['c'] = sigma['C'] = 1;
      sigma['g'] = sigma['G'] = 2;
      sigma['t'] = sigma['T'] = 3;

      std::vector<char> matrix(5 * 5, s);
      for (unsigned int i = 0; i < 5; ++i) matrix[(5 * i) + i] = m;

      return scoring_matrix(sigma, matrix);
  }; // make_dna_sm

  bool read_file_sm(const std::string& name, scoring_matrix& sub) {
      std::ifstream f(name.c_str());
      if (!f) return false;

      unsigned char sigma[256];
      std::string buf;

      // read comments
      while (!f.eof()) {
	  buf = "";
	  std::getline(f, buf);
	  if ((buf.empty() == true) || (buf[0] != '#')) break;
      } // while

      if (buf.empty() == true) return false;

      // parse column header
      std::string head = buf;
      head.erase(std::remove(head.begin(), head.end(), ' '), head.end());

      unsigned int len = head.size();
      if (head[head.size() - 1] != '*') return false;

      std::vector<char> matrix(len * len, 0);

      for (unsigned int i = 0; i < 256; ++i) sigma[i] = len - 1;
      for (unsigned int i = 0; i < len - 1; ++i) sigma[head[i]] = i;

      // read matrix
      for (unsigned int i = 0; i < len; ++i) {
	  char id = 0;
	  int val;

	  f >> id;
	  if (!f || (id != head[i])) return false;

	  for (unsigned int j = 0; j < len; ++j) {
	      f >> val;
	      if (!f) return false;
	      unsigned int pos0 = sigma[head[i]];
	      unsigned int pos1 = sigma[head[j]];
	      matrix[pos0 * len + pos1] = val;
	  } // for j
      } // for i

      f.close();

      sub = scoring_matrix(sigma, matrix);
      return true;
  } // read_scoring_matrix


  /** Class: global_alignment
   *  Functor implementing memory-efficient global pairwise sequence alignment
   *  with affine gap penalty.
   */
  class global_alignment : public sequence_compare<global_alignment> {
  public:
      /** Constructor: global_alignment
       *
       *  Parameter:
       *  m - Match score (some positive number).
       *  s - Substitution penalty (usually negative number).
       *  g - Gap opening penalty (negative number).
       *  h - Gap extension penalty (negative number).
       */
      explicit global_alignment(int m = 0, int s = 0, int g = 0, int h = 0)
	  : sub_(make_dummy_sm(m, s)), g_(g), h_(h) { }

      /** Constructor: global_alignment
       *
       *  Parameter:
       *  sm - Substitution matrix.
       *  g -  Gap opening penalty (negative number).
       *  h -  Gap extension penalty (negative number).
       */
      global_alignment(const scoring_matrix& sm, int g, int h)
	  : sub_(sm), g_(g), h_(h) { }

      /** Function: operator()
       *  Compute alignment between s0 and s1.
       *
       *  Returns:
       *  3-tuple (alignment score, alignment length without terminal gaps, number of matches).
       */
      boost::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
	  unsigned int n = s0.size() + 1;
	  unsigned int m = s1.size() + 1;

	  // S(i, j) = max{ I(i, j), D(i, j), S(i - 1, j - 1) + d(i,j) }
	  // D(i, j) = max{ D(i, j - 1), S(i, j - 1) + g } + h
	  // I(i, j) = max{ I(i - 1, j), S(i - 1, j) + g } + h

	  S_.resize(m, 0);
	  std::fill(S_.begin(), S_.end(), 0);

	  I_.resize(m, 0);
	  std::fill(I_.begin(), I_.end(), 0);

	  track_.resize(n * m);
	  std::fill(track_.begin(), track_.end(), 0);

	  for (unsigned int j = 1; j < m; ++j) {
	      track_[j] = LEFT;
	      S_[j] = I_[j] = g_ + j * h_;
	  }

	  unsigned int pos = 0;
	  int Sij = 0;

	  for (unsigned int i = 1; i < n; ++i) {
	      int Si = g_ + i * h_;
	      int Di = g_ + i * h_;

	      pos = i * m;
	      track_[pos] = TOP;

	      for (unsigned int j = 1; j < m; ++j) {
		  pos++;

		  Di = std::max(Di, Si + g_) + h_;
		  I_[j] = std::max(I_[j], S_[j] + g_) + h_;

		  Si = Sij + sub_(s0[i - 1], s1[j - 1]);

		  // default: max in Si
		  track_[pos] = DIAG;

		  if (Di < I_[j]) {
		      if (Si < I_[j]) {
			  // max in I_[j]
			  Si = I_[j];
			  track_[pos] = TOP;
		      }
		  } else {
		      if (Si < Di) {
			  // max in Di
			  Si = Di;
			  track_[pos] = LEFT;
		      }
		  } // if

		  Sij = S_[j];
		  S_[j] = Si;

	      } // for j

	      Sij = g_ + i * h_;
	  } // for i

	  // backtrack
	  unsigned int i = n - 1;
	  unsigned int j = m - 1;

	  unsigned int match = 0;
	  unsigned int length = 0;

	  bool has_gap = false;
	  unsigned int sgap = 0;
	  unsigned int egap = 0;

	  has_path_ = false;
	  path_.clear();

	  while ((i > 0) || (j > 0)) {
	      switch (track_[i * m + j]) {
		case TOP:
		    --i;
		    sgap++;
		    path_.push_back('d');
		    break;

		case LEFT:
		    --j;
		    sgap++;
		    path_.push_back('i');
		    break;

		case DIAG:
		    --i;
		    --j;

		    if (s0[i] == s1[j]) {
			match++;
			path_.push_back('m');
		    } else path_.push_back('s');

		    if (has_gap == false) {
			has_gap = true;
			egap = sgap;
		    }

		    sgap = 0;
		    break;
	      } // switch

	      length++;
	  } // while

	  return boost::make_tuple(S_.back(), length - sgap - egap, match);
      } // operator()

      /** Function: path
       *  Return the edit path of the last alignment.
       *
       *  Returns:
       *  Edit path where 'i' means insert gap in s0, 'd' is deletion,
       *  's' is substitution, and 'm' is match.
       */
      std::string path() {
	  if (has_path_ == false) {
	      std::reverse(path_.begin(), path_.end());
	      has_path_ = true;
	  }
	  return path_;
      } // path


  private:
      enum { TOP, LEFT, DIAG };

      bool has_path_;
      std::string path_;
      std::vector<unsigned char> track_;

      std::vector<int> S_;
      std::vector<int> I_;

      scoring_matrix sub_;

      int g_;
      int h_;

  }; // class global_alignment


  /** Class: cfe_global_alignment
   *  Functor implementing memory-efficient global pairwise sequence alignment
   *  with cost-free end gaps affine gap penalty. Cost-free end gaps are
   *  allowed only in one sequence.
   */
  class cfe_global_alignment : public sequence_compare<global_alignment> {
  public:
      /** Constructor: cfe_global_alignment
       *
       *  Parameter:
       *  m - Match score (some positive number).
       *  s - Substitution penalty (usually negative number).
       *  g - Gap opening penalty (negative number).
       *  h - Gap extension penalty (negative number).
       */
      explicit cfe_global_alignment(int m = 0, int s = 0, int g = 0, int h = 0)
	  : sub_(make_dummy_sm(m, s)), g_(g), h_(h) { }

      /** Constructor: cfe_global_alignment
       *
       *  Parameter:
       *  sm - Substitution matrix.
       *  g -  Gap opening penalty (negative number).
       *  h -  Gap extension penalty (negative number).
       */
      cfe_global_alignment(const scoring_matrix& sm, int g, int h)
	  : sub_(sm), g_(g), h_(h) { }

      /** Function: operator()
       *  Compute alignment between s0 and s1. s0 is assumed to be a shorter query sequence
       *  in which end-gaps are free (i.e. s0 is contained in s1).
       *
       *  Returns:
       *  3-tuple (alignment score, alignment length without terminal gaps, number of matches).
       */
      boost::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
	  unsigned int n = s0.size() + 1;
	  unsigned int m = s1.size() + 1;

	  // S(i, j) = max{ I(i, j), D(i, j), S(i - 1, j - 1) + d(i,j) }
	  // D(i, j) = max{ D(i, j - 1), S(i, j - 1) + g } + h
	  // I(i, j) = max{ I(i - 1, j), S(i - 1, j) + g } + h

	  S_.resize(m, 0);
	  std::fill(S_.begin(), S_.end(), 0);

	  I_.resize(m, 0);
	  std::fill(I_.begin(), I_.end(), 0);

	  track_.resize(n * m);
	  std::fill(track_.begin(), track_.end(), 0);

	  for (unsigned int j = 1; j < m; ++j) {
	      I_[j] = g_ + j * h_;
	      track_[j] = LEFT;
	  }

	  unsigned int pos = 0;
	  int Sij = 0;

	  for (unsigned int i = 1; i < n; ++i) {
	      int Si = g_ + i * h_;
	      int Di = g_ + i * h_;

	      pos = i * m;
	      track_[pos] = TOP;

	      for (unsigned int j = 1; j < m; ++j) {
		  pos++;

		  Di = std::max(Di, Si + g_) + h_;
		  I_[j] = std::max(I_[j], S_[j] + g_) + h_;

		  Si = Sij + sub_(s0[i - 1], s1[j - 1]);

		  // default: max in Si
		  track_[pos] = DIAG;

		  if (Di < I_[j]) {
		      if (Si < I_[j]) {
			  // max in I_[j]
			  Si = I_[j];
			  track_[pos] = TOP;
		      }
		  } else {
		      if (Si < Di) {
			  // max in Di
			  Si = Di;
			  track_[pos] = LEFT;
		      }
		  } // if

		  Sij = S_[j];
		  S_[j] = Si;

	      } // for j

	      Sij = g_ + i * h_;
	  } // for i

	  // backtrack
	  unsigned int i = n - 1;
	  unsigned int j = std::max_element(S_.begin() + 1, S_.end()) - S_.begin();

	  int score = S_[j];

	  unsigned int match = 0;
	  unsigned int length = (m - 1) - j;

	  bool has_gap = true;
	  unsigned int sgap = 0;
	  unsigned int egap = length;

	  has_path_ = false;
	  path_ = std::string(length, 'i');

	  while ((i > 0) || (j > 0)) {
	      switch (track_[i * m + j]) {
		case TOP:
		    --i;
		    sgap++;
		    path_.push_back('d');
		    break;

		case LEFT:
		    --j;
		    sgap++;
		    path_.push_back('i');
		    break;

		case DIAG:
		    --i;
		    --j;

		    if (s0[i] == s1[j]) {
			match++;
			path_.push_back('m');
		    } else path_.push_back('s');

		    if (has_gap == false) {
			has_gap = true;
			egap = sgap;
		    }

		    sgap = 0;
		    break;
	      } // switch

	      length++;
	  } // while

	  return boost::make_tuple(score, length - sgap - egap, match);
      } // operator()

      /** Function: path
       *  Return the edit path of the last alignment.
       *
       *  Returns:
       *  Edit path where 'i' means insert gap in s0, 'd' is deletion,
       *  's' is substitution, and 'm' is match.
       */
      std::string path() {
	  if (has_path_ == false) {
	      std::reverse(path_.begin(), path_.end());
	      has_path_ = true;
	  }
	  return path_;
      } // path


  private:
      enum { TOP = 0, LEFT = 1, DIAG = 2 };

      bool has_path_;
      std::string path_;
      std::vector<unsigned char> track_;

      std::vector<int> S_;
      std::vector<int> I_;

      scoring_matrix sub_;

      int g_;
      int h_;

  }; // class cfe_global_alignment


  /** Class: d2
   *  Functor to compute the d2 distance.
   *  No thread safety guarantees.
   */
  class d2 : public sequence_compare<d2> {
  public:
      /** Constructor: d2
       *
       *  Parameter:
       *  k - kmer length.
       *  isdna - assume that input sequences are DNA/RNA.
       */
      explicit d2(unsigned int k = 0, bool isdna = true) : k_(k), isdna_(isdna) { }

      /** Function: operator()
       *  Compute d2 score between s0 and s1.
       *
       *  Returns:
       *  3-tuple (d2 score, number of unique kmers in s0, number of unique kmers in s1).
       */
      boost::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
	  if ((s0.size() < k_) || (s1.size() < k_)) return boost::make_tuple(-1, -1, -1);

	  if (isdna_ == true) {
	      dC_(s0, k_, dcount0_);
	      dC_(s1, k_, dcount1_);
	      int S = detail::count_distance(dcount0_.begin(), dcount0_.end(),
					     dcount1_.begin(), dcount1_.end());
	      return boost::make_tuple(S, dcount0_.size(), dcount1_.size());
	  } else {
	      detail::general_kmer_count(s0, k_, count0_);
	      detail::general_kmer_count(s1, k_, count1_);
	      int S = detail::count_distance(count0_.begin(), count0_.end(),
					     count1_.begin(), count1_.end());
	      return boost::make_tuple(S, count0_.size(), count1_.size());
	  }

	  return boost::make_tuple(-1, -1, -1);
      } // operator()

      /** Function: operator()
       *  Compute d2 score between s0 and s1, where s0 is the sequence
       *  passed to the previous call of the binary version of this operator.
       *
       *  Returns:
       *  3-tuple (d2 score, number of unique kmers in s0, number of unique kmers in s1).
       */
      boost::tuple<int, int, int> operator()(const std::string& s1) {
	  if (s1.size() < k_) return boost::make_tuple(-1, -1, -1);

	  if (isdna_ == true) {
	      dC_(s1, k_, dcount1_);
	      int S = detail::count_distance(dcount0_.begin(), dcount0_.end(),
					     dcount1_.begin(), dcount1_.end());
	      return boost::make_tuple(S, dcount0_.size(), dcount1_.size());
	  } else {
	      detail::general_kmer_count(s1, k_, count1_);
	      int S = detail::count_distance(count0_.begin(), count0_.end(),
					     count1_.begin(), count1_.end());
	      return boost::make_tuple(S, count0_.size(), count1_.size());
	  }

	  return boost::make_tuple(-1, -1, -1);
      } // operator()

  private:
      unsigned int k_;
      bool isdna_;

      std::map<unsigned long long int, unsigned int> dcount0_;
      std::map<unsigned long long int, unsigned int> dcount1_;

      std::map<std::string, unsigned int> count0_;
      std::map<std::string, unsigned int> count1_;

      detail::dna_kmer_count dC_;

  }; // class d2


  /** Class: kmer_fraction
   *  Functor to compute the kmer fraction similarity, defined
   *  as Jaccard index between kmer spectra of sequences.
   *  No thread safety guarantees.
   */
  class kmer_fraction : public sequence_compare<kmer_fraction> {
  public:
      /** Constructor: kmer_fraction
       *
       *  Parameter:
       *  k - kmer length.
       *  isdna - assume that input sequences are DNA/RNA.
       */
      explicit kmer_fraction(unsigned int k = 0, bool isdna = true) : k_(k), isdna_(isdna) { }

      /** Function: operator()
       *  Compute kmer fraction score between s0 and s1.
       *
       *  Returns:
       *  3-tuple (score, number of kmers in s0, number of kmers in s1).
       */
      boost::tuple<int, int, int> operator()(const std::string& s0, const std::string& s1) {
	  if ((s0.size() < k_) || (s1.size() < k_)) return boost::make_tuple(-1, -1, -1);

	  if (isdna_ == true) {
	      dI_(s0, k_, dindex0_);
	      dI_(s1, k_, dindex1_);
	      int S = detail::intersection_size(dindex0_.begin(), dindex0_.end(),
						dindex1_.begin(), dindex1_.end(),
						std::less<unsigned long long int>());
	      return boost::make_tuple(S, dindex0_.size(), dindex1_.size());
	  } else {
	      detail::general_kmer_index(s0, k_, index0_);
	      detail::general_kmer_index(s1, k_, index1_);
	      int S = detail::intersection_size(index0_.begin(), index0_.end(),
						index1_.begin(), index1_.end(),
						std::less<std::string>());
	      return boost::make_tuple(S, index0_.size(), index1_.size());
	  }

	  return boost::make_tuple(-1, -1, -1);
      } // operator()

      /** Function: operator()
       *  Compute kmer fraction score between s0 and s1, where s0 is the sequence
       *  passed to the previous call of the binary version of this operator.
       *
       *  Returns:
       *  3-tuple (score, number of kmers in s0, number of kmers in s1).
       */
      boost::tuple<int, int, int> operator()(const std::string& s1) {
	  if (s1.size() < k_) return boost::make_tuple(-1, -1, -1);

	  if (isdna_ == true) {
	      dI_(s1, k_, dindex1_);
	      int S = detail::intersection_size(dindex0_.begin(), dindex0_.end(),
						dindex1_.begin(), dindex1_.end(),
						std::less<unsigned long long int>());
	      return boost::make_tuple(S, dindex0_.size(), dindex1_.size());
	  } else {
	      detail::general_kmer_index(s1, k_, index1_);
	      int S = detail::intersection_size(index0_.begin(), index0_.end(),
						index1_.begin(), index1_.end(),
						std::less<std::string>());
	      return boost::make_tuple(S, index0_.size(), index1_.size());
	  }

	  return boost::make_tuple(-1, -1, -1);
      } // operator()

  private:
      unsigned int k_;
      bool isdna_;

      std::vector<unsigned long long int> dindex0_;
      std::vector<unsigned long long int> dindex1_;

      std::vector<std::string> index0_;
      std::vector<std::string> index1_;

      detail::dna_kmer_index dI_;

  }; // class kmer_fraction

} // namespace bio

#endif // SEQUENCE_COMPARE_HPP
