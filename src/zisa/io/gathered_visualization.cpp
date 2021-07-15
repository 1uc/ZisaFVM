// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/io/gathered_visualization.hpp>

#include <zisa/math/permutation.hpp>

namespace zisa {

GatheredVisualization::GatheredVisualization(
    std::unique_ptr<AllVariablesGatherer> all_vars_gatherer,
    std::shared_ptr<Permutation> permutation,
    std::shared_ptr<Visualization> visualization,
    const AllVariablesDimensions &all_var_dims)
    : gatherer(std::move(all_vars_gatherer)),
      permutation(std::move(permutation)),
      visualization(std::move(visualization)) {

  LOG_ERR_IF(this->gatherer == nullptr, "Received a `nullptr`.");

  if (this->gatherer->is_this_rank_gathering()) {
    LOG_ERR_IF(this->permutation == nullptr, "Received a `nullptr`.");
    LOG_ERR_IF(this->visualization == nullptr, "Received a `nullptr`.");

    buffer = AllVariables(all_var_dims);
  }
}

GatheredVisualization::~GatheredVisualization() {
  GatheredVisualization::do_wait();
}

template <class Vis>
void GatheredVisualization::gather_and_visualize(
    const AllVariables &all_variables, const Vis &vis) {

  if (gatherer->is_this_rank_gathering()) {
    wait();

    gatherer->copy_local_patch(buffer, all_variables);

    job = std::make_unique<std::thread>([this, vis]() {
      gatherer->receive(buffer);

      apply_permutation(array_view(buffer.cvars), *permutation);
      apply_permutation(array_view(buffer.avars), *permutation);

      vis(buffer);
    });
  } else {
    gatherer->send(all_variables);
  }
}

void GatheredVisualization::do_visualization(
    const AllVariables &all_variables,
    const SimulationClock &simulation_clock) {

  auto vis = [this, &simulation_clock](const AllVariables &full_vars) {
    (*visualization)(full_vars, simulation_clock);
  };

  gather_and_visualize(all_variables, vis);
}

void GatheredVisualization::do_steady_state(const AllVariables &all_variables) {
  auto vis = [this](const AllVariables &full_vars) {
    visualization->steady_state(full_vars);
  };

  gather_and_visualize(all_variables, vis);
}

void GatheredVisualization::do_wait() {
  if (job != nullptr && job->joinable()) {
    job->join();
  }
}

}