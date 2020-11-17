#include <zisa/mpi/io/gathered_visualization_factory.hpp>

namespace zisa {
std::shared_ptr<GatheredVisualization> make_gathered_visualization(
    std::shared_ptr<GatheredVisInfo> vis_info,
    std::shared_ptr<MPISingleNodeArrayGathererFactory> gatherer_factory,
    std::shared_ptr<Visualization> visualization,
    AllVariablesDimensions all_var_dims) {

  auto n_local_cells = vis_info->n_local_cells;
  auto n_vis_cells = vis_info->n_vis_cells();

  auto cvars_gatherer = gatherer_factory->template create_pointer<double, 2>();
  auto avars_gatherer = gatherer_factory->template create_pointer<double, 2>();

  auto all_var_gatherer
      = std::make_unique<AllVariablesGatherer>(std::move(cvars_gatherer),
                                               std::move(avars_gatherer),
                                               n_local_cells,
                                               n_vis_cells);

  if (vis_info->h5_comm != MPI_COMM_NULL) {
    all_var_dims.n_cells = n_vis_cells;
    return std::make_shared<GatheredVisualization>(std::move(all_var_gatherer),
                                                   vis_info->permutation,
                                                   visualization,
                                                   all_var_dims);
  } else {
    return std::make_shared<GatheredVisualization>(
        std::move(all_var_gatherer),
        nullptr,
        nullptr,
        AllVariablesDimensions{0, 0, 0});
  }
}
}