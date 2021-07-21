// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/mpi/parallelization/mpi_single_node_array_gatherer.hpp>

#include <zisa/mpi/mpi_tag_constants.hpp>

namespace zisa {

std::shared_ptr<MPISingleNodeArrayGathererFactory>
make_mpi_single_node_array_gatherer_factory(const GatheredVisInfo &vis_info) {
  auto vis_comm = vis_info.vis_comm;
  auto vis_tag = ZISA_MPI_TAG_GATHERED_VIS;

  auto h5_comm = vis_info.h5_comm;

  if (h5_comm != MPI_COMM_NULL) {
    auto darray_info = make_distributed_array_info(vis_info.vis_boundaries);
    return std::make_shared<MPISingleNodeArrayGathererFactory>(
        darray_info, vis_comm, vis_tag + 1);
  } else {
    return std::make_shared<MPISingleNodeArrayGathererFactory>(
        nullptr, vis_comm, vis_tag + 1);
  }
}

}
