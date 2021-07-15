// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/io/load_snapshot.hpp>

namespace zisa {

LoadSnapshot::LoadSnapshot(std::shared_ptr<FNG> fng) : fng(std::move(fng)) {}

void LoadSnapshot::do_load(AllVariables &all_variables,
                           SimulationClock &simulation_clock) {
  auto reader = pick_reader(fng->next_name());
  auto [time, n_steps] = load_full_state(*reader, all_variables);
  simulation_clock.advance_to(time, n_steps);
}

std::unique_ptr<HDF5Reader>
SerialLoadSnapshot::pick_reader(const std::string &file_name) {
  return std::make_unique<HDF5SerialReader>(file_name);
}

}