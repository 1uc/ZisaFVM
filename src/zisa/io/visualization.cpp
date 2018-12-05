/* Implementation for visualizing the flow.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-09-29
 */
#include <zisa/io/visualization.hpp>

namespace zisa {

void Visualization::operator()(const AllVariables &all_variables,
                               const SimulationClock &simulation_clock) {
  do_visualization(all_variables, simulation_clock);
}

} // namespace zisa
