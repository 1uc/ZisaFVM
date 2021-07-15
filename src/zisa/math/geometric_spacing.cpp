// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/basic_functions.hpp>
#include <zisa/math/geometric_spacing.hpp>

namespace zisa {
array<double, 1>
geometric_spacing(double x_end, double growth_factor, int_t n_points) {

  // L = dx * \sum_k a^k = dx * (1 - a^n) / (1 - a)

  double a = growth_factor;
  double dx = (1.0 - a) / (1.0 - zisa::pow(a, double(n_points))) * x_end;

  auto spacing = array<double, 1>(n_points);
  spacing[0] = 0.0;

  double ak = 1.0;
  for (int_t k = 1; k < n_points; ++k) {
    ak *= a;
    spacing[k] = spacing[k - 1] + dx * ak;
  }

  spacing[n_points - 1] = x_end;

  return spacing;
}
}