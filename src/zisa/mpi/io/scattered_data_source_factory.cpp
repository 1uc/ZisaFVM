// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/mpi/io/scattered_data_source_factory.hpp>

namespace zisa {

std::shared_ptr<ScatteredDataSource> make_scattered_data_source(
    std::shared_ptr<GatheredVisInfo> vis_info,
    std::shared_ptr<MPISingleNodeArrayScattererFactory> scatterer_factory,
    std::shared_ptr<DataSource> data_source,
    std::shared_ptr<HaloExchange> halo_exchange,
    AllVariablesDimensions all_var_dims) {

  auto n_local_cells = vis_info->n_local_cells;
  auto n_vis_cells = vis_info->n_vis_cells();

  auto cvars_scatterer
      = scatterer_factory->template create_pointer<double, 2>();
  auto avars_scatterer
      = scatterer_factory->template create_pointer<double, 2>();

  auto all_var_scatterer
      = std::make_unique<AllVariablesScatterer>(std::move(cvars_scatterer),
                                                std::move(avars_scatterer),
                                                n_local_cells,
                                                n_vis_cells);

  if (vis_info->h5_comm != MPI_COMM_NULL) {
    all_var_dims.n_cells = n_vis_cells;
    return std::make_shared<ScatteredDataSource>(std::move(all_var_scatterer),
                                                 vis_info->permutation,
                                                 std::move(data_source),
                                                 std::move(halo_exchange),
                                                 all_var_dims);
  } else {
    return std::make_shared<ScatteredDataSource>(
        std::move(all_var_scatterer),
        nullptr,
        nullptr,
        std::move(halo_exchange),
        AllVariablesDimensions{0, 0, 0});
  }
}

}
