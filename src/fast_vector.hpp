/***
 *  $Id$
 **
 *  File: fast_vector.hpp
 *  Created: Jan 07, 2014
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2014 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef FAST_VECTOR_HPP
#define FAST_VECTOR_HPP

#include <cstdlib>
#include <cstring>
#include <iterator>
#include <new>

template <typename T> class fast_vector {
public:
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef std::size_t size_type;
    typedef std::size_t difference_type;
    typedef T*  pointer;
    typedef const T* const_pointer;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    explicit fast_vector(size_type n = 0) : sz_(0), n_(0), buf_(0) {
        if (n > 0) reserve(n);
        n_ = n;
    } // fast_vector

    fast_vector(const fast_vector<value_type>& v) : sz_(0), n_(0), buf_(0) {
        reserve(v.n_);
        n_ = v.n_;
        std::memcpy(buf_, v.buf_, n_ * sizeof(value_type));
    } // fast_vector

    ~fast_vector() {
        if (buf_ != 0) std::free(buf_);
        sz_ = 0;
        n_ = 0;
        buf_ = 0;
    } // ~fast_vector

    fast_vector<value_type>& operator=(const fast_vector<value_type>& v) {
        if (this == &v) return *this;
        reserve(v.n_);
        n_ = v.n_;
        std::memcpy(buf_, v.buf_, n_ * sizeof(value_type));
        return *this;
    } // operator=

    iterator begin() { return buf_; }

    iterator end() { return buf_ + n_; }

    const_iterator begin() const { return buf_; }

    const_iterator end() const { return buf_ + n_; }

    reverse_iterator rbegin() { return reverse_iterator(begin()); }

    reverse_iterator rend() { return reverse_iterator(end()); }

    const_reverse_iterator rbegin() const { return reverse_iterator(begin()); }

    const_reverse_iterator rend() const { return reverse_iterator(end()); }

    size_type size() const { return n_; }

    size_type capacity() const { return sz_; }

    bool empty() const { return (n_ == 0); }

    reference front() { return buf_[0]; }

    const_reference front() const { return buf_[0]; }

    reference back() { return buf_[n_ - 1]; }

    const_reference back() const { return buf_[n_ - 1]; }

    reference operator[](size_type n) { return buf_[n]; }

    const_reference operator[](size_type n) const { return buf_[n]; }

    void clear() { resize(0); }

    void reserve(size_type n) {
        if (sz_ < n) {
            sz_ = n;
            buf_ = static_cast<pointer>(std::realloc(buf_, sz_ * sizeof(value_type)));
            if (buf_ == 0) throw std::bad_alloc();
        }
    } // reserve

    void resize(size_type n) {
        sz_ = n;
        n_ = n;
        buf_ = static_cast<pointer>(std::realloc(buf_, sz_ * sizeof(value_type)));
        if (buf_ == 0) throw std::bad_alloc();
    } // resize

    void push_back(const value_type& t) {
        if (sz_ <= n_) {
            sz_ += BLOCK;
            buf_ = static_cast<pointer>(std::realloc(buf_, sz_ * sizeof(value_type)));
        }
        if (buf_ == 0) throw std::bad_alloc();
        buf_[n_] = t;
        n_++;
    } // push_back

    void swap(fast_vector<value_type>& v) {
        std::swap(sz_, v.sz_);
        std::swap(n_, v.n_);
        std::swap(buf_, v.buf_);
    } // swap

private:
    enum { BLOCK = 4096 };
    size_type sz_;
    size_type n_;
    pointer buf_;

}; // class fast_vector

#endif // FAST_VECTOR_HPP
