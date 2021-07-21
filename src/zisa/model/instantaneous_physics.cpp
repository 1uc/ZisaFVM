// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/model/instantaneous_physics.hpp>

namespace zisa {
void NoInstantaneousPhysics::compute(const zisa::SimulationClock &,
                                     zisa::AllVariables &) {
  return;
}
}