/***
 *  $Id$
 **
 *  File: sequence_compact.hpp
 *  Created: Jul 10, 2015
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@hush.com>
 *  Copyright (c) 2015 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef SEQUENCE_COMPACT_HPP
#define SEQUENCE_COMPACT_HPP

#include <algorithm>
#include <cmath>

#include <jaz/algorithm.hpp>


template <typename Hash> class to_move {
public:
    to_move(Hash hash, int root, int size) : hash_(hash), root_(root), size_(size) { }

    template <typename T>
    bool operator()(const T& t) const { return (hash_(t) % size_ != root_); }

private:
    Hash hash_;
    int root_;
    int size_;

}; // to_move


template <typename Sequence, typename Hash, typename LTComp, typename Oper, typename EqComp>
void gather_compact(Sequence& seq_in, Sequence& seq_out,
                    Hash hash, LTComp ltcomp, Oper op, EqComp eqcomp,
                    MPI_Datatype Type, int root, MPI_Comm Comm) {
    typedef typename Sequence::value_type value_type;
    typedef typename Sequence::iterator iterator;

    int size = 0;
    int rank = 0;

    MPI_Comm_size(Comm, &size);
    MPI_Comm_rank(Comm, &rank);

    int p = size;
    if (p == 1) return;

    // analyze local buffer to select objects for root
    iterator iter = std::stable_partition(seq_in.begin(), seq_in.end(), to_move<Hash>(hash, root, size));
    int pos = iter - seq_in.begin();
    int m = (seq_in.end() - iter);

    // gather and compact
    std::vector<int> sz(size);
    std::vector<int> displ;

    int S = 0;

    MPI_Gather(&m, 1, MPI_INT, &sz[0], 1, MPI_INT, root, Comm);

    if (rank == root) {
        S = std::accumulate(sz.begin(), sz.end(), 0);
        seq_out.resize(S);

        displ.resize(size, 0);
        std::partial_sum(sz.begin(), sz.end() - 1, displ.begin() + 1);

        MPI_Gatherv(&seq_in[pos], m, Type, &seq_out[0], &sz[0], &displ[0], Type, root, Comm);

        std::sort(seq_out.begin(), seq_out.end(), ltcomp);
        seq_out.resize(jaz::compact(seq_out.begin(), seq_out.end(), op, eqcomp) - seq_out.begin());
    } else {
        MPI_Gatherv(&seq_in[pos], m, Type, &seq_out[0], &sz[0], &displ[0], Type, root, Comm);
    }

    seq_in.resize(pos);
} // gather_compact

template <typename Sequence, typename Hash, typename Oper>
void gather_compact(Sequence& seq_in, Sequence& seq_out, Hash hash, Oper op,
                    MPI_Datatype Type, int root, MPI_Comm Comm) {
    typedef typename Sequence::value_type value_type;
    gather_compact(seq_in, seq_out, hash, std::less<value_type>(), op, std::equal_to<value_type>(),
                   Type, root, Comm);
} // gather_compact


template <typename Sequence, typename Hash, typename LTComp, typename Oper, typename EqComp>
void tree_compact(Sequence& seq_in, Sequence& seq_out,
                  Hash hash, LTComp ltcomp, Oper op, EqComp eqcomp,
                  MPI_Datatype Type, int root, MPI_Comm Comm) {
    typedef typename Sequence::value_type value_type;
    typedef typename Sequence::iterator iterator;

    int size = 0;
    int rank = 0;

    MPI_Comm_size(Comm, &size);
    MPI_Comm_rank(Comm, &rank);

    int p = size;
    if (p == 1) return;

    // analyze local buffer to select objects for root
    iterator iter = std::stable_partition(seq_in.begin(), seq_in.end(), to_move<Hash>(hash, root, size));
    int m = (seq_in.end() - iter);

    seq_out.resize(m);
    std::copy(iter, seq_in.end(), seq_out.begin());

    seq_in.resize(seq_in.size() - m);

    // execute reduction over the virtual tree
    const int TAG = 333;
    MPI_Status stat;

    int r = ceil(log2(size));
    int d = 1;

    int id = rank;

    if (rank == 0) id = root;
    if (rank == root) id = 0;

    for (int i = 0; i < r; ++i, d = d << 1) {
        if (id % d == 0) {
            if (id % (d << 1) == 0) {
                int src = id + d;

                if (src == root) src = 0;
                else if (src == 0) src = root;

                if (src < size) {
                    int sz = 0;
                    MPI_Recv(&sz, 1, MPI_INT, src, TAG, Comm, &stat);

                    int pos = seq_out.size();
                    seq_out.resize(pos + sz);

                    MPI_Recv(&seq_out[pos], sz, Type, src, TAG, Comm, &stat);

                    // merging
                    std::inplace_merge(seq_out.begin(), seq_out.begin() + pos, seq_out.end(), ltcomp);
                    seq_out.resize(jaz::compact(seq_out.begin(), seq_out.end(), op, eqcomp) - seq_out.begin());
                }
            } else {
                int dest = id - d;

                if (dest == root) dest = 0;
                else if (dest == 0) dest = root;

                int sz = seq_out.size();

                MPI_Send(&sz, 1, MPI_INT, dest, TAG, Comm);
                MPI_Send(&seq_out[0], sz, Type, dest, TAG, Comm);

                Sequence().swap(seq_out);
            }
        } // if id
    } // for i
} // tree_compact

template <typename Sequence, typename Hash, typename Oper>
void tree_compact(Sequence& seq_in, Sequence& seq_out, Hash hash, Oper op,
                      MPI_Datatype Type, int root, MPI_Comm Comm) {
    typedef typename Sequence::value_type value_type;
    tree_compact(seq_in, seq_out, hash, std::less<value_type>(), op, std::equal_to<value_type>(),
                 Type, root, Comm);
} // tree_compact

#endif // SEQUENCE_COMPACT_HPP
