// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_DECL_HPP_XCJQPO
#define ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_DECL_HPP_XCJQPO

#include <zisa/config.hpp>

#include <zisa/mpi/io/gathered_vis_info.hpp>
#include <zisa/mpi/mpi.hpp>
#include <zisa/mpi/parallelization/distributed_array_info.hpp>
#include <zisa/parallelization/array_scatterer.hpp>

namespace zisa {

template <class T, int n_dims>
class MPISingleNodeArrayScatterer : public ArrayScatterer<T, n_dims> {
private:
  using super = ArrayScatterer<T, n_dims>;

protected:
  using view_t = typename super::view_t;
  using const_view_t = typename super::const_view_t;

public:
  MPISingleNodeArrayScatterer(std::shared_ptr<DistributedArrayInfo> array_info,
                              MPI_Comm mpi_comm,
                              int mpi_tag);

  void send(const const_view_t &const_view) const override;
  void receive(const view_t &view) const override;
  void copy_local_patch(const view_t &local,
                        const const_view_t &global) const override;

  bool is_this_rank_scattering() const override;

private:
  std::shared_ptr<DistributedArrayInfo> array_info;

  int tag;
  int rank;
  int comm_size;
  int scatter_rank = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
};

class MPISingleNodeArrayScattererFactory {
public:
  MPISingleNodeArrayScattererFactory(
      std::shared_ptr<DistributedArrayInfo> array_info,
      MPI_Comm comm,
      int mpi_base_tag);

  template <class T, int n_dims>
  MPISingleNodeArrayScatterer<T, n_dims> create_object() {
    return MPISingleNodeArrayScatterer<T, n_dims>(
        array_info, comm, current_mpi_tag++);
  }

  template <class T, int n_dims>
  std::unique_ptr<MPISingleNodeArrayScatterer<T, n_dims>> create_pointer() {
    return std::make_unique<MPISingleNodeArrayScatterer<T, n_dims>>(
        array_info, comm, current_mpi_tag++);
  }

private:
  std::shared_ptr<DistributedArrayInfo> array_info;
  MPI_Comm comm;
  int current_mpi_tag;
};

std::shared_ptr<MPISingleNodeArrayScattererFactory>
make_mpi_single_node_array_scatterer_factory(const GatheredVisInfo &vis_info);

}

#endif // ZISA_MPI_SINGLE_NODE_ARRAY_SCATTERER_DECL_HPP
