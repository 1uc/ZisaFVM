// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/io/no_visualization.hpp>

namespace zisa {

void NoVisualization::do_visualization(const AllVariables &,
                                       const SimulationClock &) {
  return;
}

} // namespace zisa
