#include <zisa/parallelization/mpi_all_variables_gatherer.hpp>

namespace zisa {

std::unique_ptr<AllVariablesGatherer> make_mpi_all_variables_gatherer(
    int base_tag, const MPI_Comm &comm, const array<int_t, 1> &boundaries) {

  auto array_info = std::make_shared<DistributedArrayInfo>(boundaries);
  auto cvars_gatherer = std::make_unique<MPISingleNodeArrayGatherer<double, 2>>(
      array_info, comm, base_tag);

  auto avars_gatherer = std::make_unique<MPISingleNodeArrayGatherer<double, 2>>(
      array_info, comm, base_tag + 1);

  auto rank = zisa::mpi::rank(comm);
  auto n_owned_cells = boundaries[rank + 1] - boundaries[rank];
  auto n_cells_total = boundaries[boundaries.size() - 1] - boundaries[0];

  return std::make_unique<AllVariablesGatherer>(std::move(cvars_gatherer),
                                                std::move(avars_gatherer),
                                                n_owned_cells,
                                                n_cells_total);
}

}
