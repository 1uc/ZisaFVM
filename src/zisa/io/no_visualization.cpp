#include <zisa/io/no_visualization.hpp>

namespace zisa {

void NoVisualization::do_visualization(const AllVariables &,
                                       const SimulationClock &) {
  return;
}

} // namespace zisa
