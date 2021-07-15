// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/mpi/io/gathered_vis_info.hpp>

#include <algorithm>
#include <numeric>
#include <zisa/math/comparison.hpp>
#include <zisa/mpi/mpi_tag_constants.hpp>

namespace zisa {

int_t GatheredVisInfo::n_vis_cells() const {
  if (h5_comm != MPI_COMM_NULL) {
    return vis_boundaries[vis_boundaries.size() - 1];
  } else {
    return 0;
  }
}

std::shared_ptr<GatheredVisInfo> make_gathered_vis_info(
    MPI_Comm world_comm, const DistributedGrid &dgrid, int n_writers) {

  auto world_rank = zisa::mpi::rank(world_comm);
  auto world_size = zisa::mpi::size(world_comm);

  auto n_local_cells = int_t(std::count(
      dgrid.partition.begin(), dgrid.partition.end(), int_t(world_rank)));

  auto n_vis_pes = zisa::min(n_writers, world_size);
  auto pes_per_block = (world_size + n_vis_pes - 1) / n_vis_pes;
  auto vis_id = world_rank / pes_per_block;

  auto vis_comm = zisa::mpi::comm_split(world_comm, vis_id, world_rank);
  auto vis_size = zisa::mpi::size(vis_comm);
  auto vis_rank = zisa::mpi::rank(vis_comm);
  auto vis_tag = ZISA_MPI_TAG_GATHERED_VIS;

  auto h5_comm = zisa::mpi::comm_split(
      world_comm, vis_rank == 0 ? 0 : MPI_UNDEFINED, world_rank);

  {
    std::string name = "VisComm-";
    name += std::to_string(vis_id);
    MPI_Comm_set_name(vis_comm, name.c_str());
  }

  if (h5_comm != MPI_COMM_NULL) {
    auto vis_cells_per_pe = array<int_t, 1>(shape_t<1>{vis_size});
    vis_cells_per_pe[0] = n_local_cells;
    zisa::mpi::gather(array_view(vis_cells_per_pe), 0, vis_comm);
    auto n_vis_cells = std::accumulate(
        vis_cells_per_pe.begin(), vis_cells_per_pe.end(), int_t(0));

    auto gids = array<int_t, 1>(shape_t<1>{n_vis_cells});
    auto vis_boundaries = array<int_t, 1>(shape_t<1>{vis_size + 1});
    vis_boundaries[0] = 0;
    for (int i = 0; i < vis_size; ++i) {
      vis_boundaries[i + 1] = vis_boundaries[i] + vis_cells_per_pe[i];
    }

    for (int_t i = 0; i < n_local_cells; ++i) {
      gids[i] = dgrid.global_cell_indices[i];
    }

    for (int p = 1; p < vis_size; ++p) {
      auto i0 = vis_boundaries[p];
      auto i1 = vis_boundaries[p + 1];
      zisa::mpi::recv(slice(array_view(gids), i0, i1), p, vis_tag, vis_comm);
    }

    auto sigma_array = array<int_t, 1>(gids.shape());
    for (int_t i = 0; i < gids.shape(0); ++i) {
      sigma_array[i] = i;
    }
    std::sort(sigma_array.begin(),
              sigma_array.end(),
              [&gids](size_t i, size_t j) { return gids[i] < gids[j]; });

    auto sigma = std::make_shared<Permutation>(
        factor_permutation(array_const_view(sigma_array)));
    apply_permutation(array_view(gids), *sigma);

    return std::make_shared<GatheredVisInfo>(
        GatheredVisInfo{std::move(vis_boundaries),
                        std::move(gids),
                        std::move(sigma),
                        n_local_cells,
                        vis_comm,
                        h5_comm});

  } else {
    zisa::mpi::gather(array_view(shape_t<1>{1}, &n_local_cells), 0, vis_comm);
    zisa::mpi::send(const_slice(array_const_view(dgrid.global_cell_indices),
                                0,
                                n_local_cells),
                    0,
                    vis_tag,
                    vis_comm);

    return std::make_shared<GatheredVisInfo>(GatheredVisInfo{array<int_t, 1>(),
                                                             array<int_t, 1>(),
                                                             nullptr,
                                                             n_local_cells,
                                                             vis_comm,
                                                             h5_comm});
  }
}

std::shared_ptr<HDF5UnstructuredFileDimensions>
make_hdf5_unstructured_file_dimensions(const GatheredVisInfo &vis_info) {

  if (vis_info.h5_comm != MPI_COMM_NULL) {
    const auto &gids = vis_info.vis_file_ids;
    std::vector<hsize_t> hids(gids.size());
    for (int_t i = 0; i < hids.size(); ++i) {
      hids[i] = integer_cast<hsize_t>(gids[i]);
    }

    return make_hdf5_unstructured_file_dimensions(
        vis_info.n_vis_cells(), hids, vis_info.h5_comm);
  } else {
    return std::make_shared<HDF5UnstructuredFileDimensions>(
        0, 0, std::vector<hsize_t>{}, MPI_COMM_NULL);
  }
}
}
