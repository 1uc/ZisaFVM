// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef EULER_PLOTS_H_G38NN
#define EULER_PLOTS_H_G38NN

#include <zisa/config.hpp>
#include <zisa/io/scalar_plot.hpp>
#include <zisa/io/visualization.hpp>

#include <chrono>
#include <thread>

namespace zisa {

class EulerPlots : public Visualization {
public:
  EulerPlots(const Grid &grid, const std::chrono::milliseconds &delay)
      : rho_plot(grid.vertices, grid.vertex_indices, "Density"), delay(delay) {}

protected:
  virtual void do_visualization(const AllVariables &all_variables,
                                const SimulationClock &) override {
    rho_plot([&all_variables](int_t i) { return all_variables.cvars(i, 0); });
    std::this_thread::sleep_for(delay);
  }

private:
  ScalarPlot rho_plot;
  std::chrono::milliseconds delay;
};

} // namespace zisa

#endif
