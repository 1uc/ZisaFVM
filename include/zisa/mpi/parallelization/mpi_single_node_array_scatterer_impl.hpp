// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_IMPL_HPP_CJUHWE
#define ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_IMPL_HPP_CJUHWE

#include "mpi_single_node_array_scatterer_decl.hpp"

namespace zisa {

template <class T, int n_dims>
MPISingleNodeArrayScatterer<T, n_dims>::MPISingleNodeArrayScatterer(
    std::shared_ptr<DistributedArrayInfo> array_info,
    MPI_Comm mpi_comm,
    int mpi_tag)
    : array_info(std::move(array_info)), tag(mpi_tag), comm(mpi_comm) {

  rank = zisa::mpi::rank(comm);
  comm_size = zisa::mpi::size(comm);
}

template <class T, int n_dims>
void MPISingleNodeArrayScatterer<T, n_dims>::send(
    const const_view_t &const_view) const {

  assert(rank == scatter_rank);

  std::vector<zisa::mpi::Request> requests;
  requests.reserve(integer_cast<size_t>(comm_size));

  for (int receiver_rank = 0; receiver_rank < comm_size; ++receiver_rank) {
    if (receiver_rank != rank) {
      auto i0 = (*array_info).partition[receiver_rank];
      auto i1 = (*array_info).partition[receiver_rank + 1];
      auto sub_array = const_slice(const_view, i0, i1);

      auto request = zisa::mpi::isend(sub_array, receiver_rank, tag, comm);
      requests.push_back(std::move(request));
    }
  }

  zisa::mpi::wait_all(requests);
}

template <class T, int n_dims>
void MPISingleNodeArrayScatterer<T, n_dims>::receive(const view_t &view) const {
  assert(rank != scatter_rank);
  zisa::mpi::recv(view, scatter_rank, tag, comm);
}

template <class T, int n_dims>
void MPISingleNodeArrayScatterer<T, n_dims>::copy_local_patch(
    const view_t &local, const const_view_t &global) const {
  int_t i0 = (*array_info).partition[rank];
  int_t i1 = (*array_info).partition[rank + 1];

  local.copy_data(const_slice(global, i0, i1));
}

template <class T, int n_dims>
bool MPISingleNodeArrayScatterer<T, n_dims>::is_this_rank_scattering() const {
  return rank == scatter_rank;
}

}

#endif // ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_IMPL_HPP
