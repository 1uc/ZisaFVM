// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/io/data_source.hpp>

namespace zisa {

void DataSource::operator()(AllVariables &all_variables,
                            SimulationClock &simulation_clock) {
  this->do_load(all_variables, simulation_clock);
}

}