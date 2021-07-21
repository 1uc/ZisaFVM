// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/experiments/smooth_bubble.hpp>

namespace zisa {

std::pair<std::shared_ptr<AllVariables>, std::shared_ptr<AllVariables>>
SmoothBubble::compute_initial_conditions() {
  auto all_variables
      = std::make_shared<AllVariables>(choose_all_variable_dims());
  auto &u0 = all_variables->cvars;
  auto grid = choose_grid();

  auto x_ref = XYZ{0.0, 0.0, 0.0};

  double E = 1.0;

  const auto &ic = params["experiment"]["initial_conditions"];
  double A = ic["amplitude"];
  double sigma = ic["width"];

  for (const auto &[i, tri] : triangles(*grid)) {
    double r = norm(barycenter(tri) - x_ref);
    double dE = A * zisa::exp(-zisa::pow<2>(r / sigma));
    u0(i) = cvars_t{0.1, 0.0, 0.0, 0.0, E + dE};
  }

  return {all_variables, nullptr};
}

} // namespace zisa
