// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_IMPL_HPP
#define ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_IMPL_HPP

#include "zisa/mpi/parallelization/mpi_single_node_array_gatherer_decl.hpp"

#include <zisa/utils/integer_cast.hpp>

namespace zisa {

template <class T, int n_dims>
MPISingleNodeArrayGatherer<T, n_dims>::MPISingleNodeArrayGatherer(
    std::shared_ptr<DistributedArrayInfo> array_info,
    MPI_Comm mpi_comm,
    int mpi_tag)
    : array_info(std::move(array_info)), tag(mpi_tag), comm(mpi_comm) {

  rank = zisa::mpi::rank(comm);
  comm_size = zisa::mpi::size(comm);
}

template <class T, int n_dims>
void MPISingleNodeArrayGatherer<T, n_dims>::send(
    const MPISingleNodeArrayGatherer::const_view_t &const_view) const {
  assert(rank != gather_rank);

  zisa::mpi::send(const_view, gather_rank, tag, comm);
}

template <class T, int n_dims>
void MPISingleNodeArrayGatherer<T, n_dims>::receive(
    const MPISingleNodeArrayGatherer::view_t &view) const {
  std::vector<zisa::mpi::Request> requests;
  requests.reserve(integer_cast<size_t>(comm_size));

  for (int sender_rank = 0; sender_rank < comm_size; ++sender_rank) {
    if (sender_rank != rank) {
      auto i0 = (*array_info).partition[sender_rank];
      auto i1 = (*array_info).partition[sender_rank + 1];
      auto sub_array = slice(view, i0, i1);

      auto request = zisa::mpi::irecv(sub_array, sender_rank, tag, comm);
      requests.push_back(std::move(request));
    }
  }

  zisa::mpi::wait_all(requests);
}

template <class T, int n_dims>
void MPISingleNodeArrayGatherer<T, n_dims>::copy_local_patch(
    const MPISingleNodeArrayGatherer::view_t &global,
    const MPISingleNodeArrayGatherer::const_view_t &local) const {
  int_t i0 = (*array_info).partition[rank];
  int_t i1 = (*array_info).partition[rank + 1];

  slice(global, i0, i1).copy_data(local);
}

template <class T, int n_dims>
bool MPISingleNodeArrayGatherer<T, n_dims>::is_this_rank_gathering() const {
  return rank == 0;
}
}

#endif // ZISA_MPI_SINGLE_NODE_ARRAY_GATHERER_IMPL_HPP
