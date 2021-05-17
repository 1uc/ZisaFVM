#include <zisa/io/visualization.hpp>

namespace zisa {

void Visualization::operator()(const AllVariables &all_variables,
                               const SimulationClock &simulation_clock) {
  do_visualization(all_variables, simulation_clock);
}

void Visualization::steady_state(const AllVariables &all_variables) {
  do_steady_state(all_variables);
}

void Visualization::do_steady_state(const AllVariables &) {
  LOG_ERR("You need to implement this in the subclass.");
}
void Visualization::wait() { do_wait(); }
void Visualization::do_wait() { return; }

} // namespace zisa
