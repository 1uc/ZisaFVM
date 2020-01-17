#include <zisa/io/gathered_visualization.hpp>

#include <zisa/parallelization/mpi.hpp>

namespace zisa {

GatheredVisualization::GatheredVisualization(
    std::shared_ptr<AllVariablesGatherer> all_vars_gatherer,
    std::shared_ptr<Visualization> visualization,
    const AllVariablesDimensions &all_var_dims)
    : gatherer(std::move(all_vars_gatherer)),
      visualization(std::move(visualization)) {

  buffer = AllVariables(all_var_dims);
}

GatheredVisualization::~GatheredVisualization() {
  if (job != nullptr && job->joinable()) {
    job->join();
  }
}

void GatheredVisualization::do_visualization(
    const AllVariables &all_variables,
    const SimulationClock &simulation_clock) {
  LOG_ERR_IF(all_variables.avars.shape(1) != 0,
             "Implement advected variables first.");

  if (gatherer->is_this_rank_gathering()) {
    gatherer->copy_local_patch(buffer, all_variables);

    if (job != nullptr) {
      job->join();
    }

    job = std::make_unique<std::thread>([this, &simulation_clock]() {
      gatherer->receive(buffer);
      (*visualization)(buffer, simulation_clock);
    });
  } else {
    gatherer->send(all_variables);
  }
}
}