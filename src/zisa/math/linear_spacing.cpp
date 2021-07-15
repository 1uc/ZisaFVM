// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/linear_spacing.hpp>

namespace zisa {

array<double, 1> linear_spacing(double x0, double x1, int_t n_points) {
  auto spacing = array<double, 1>(n_points);

  double dx = (x1 - x0) / double(n_points - 1);
  for (int_t i = 0; i < n_points; ++i) {
    spacing[i] = x0 + double(i) * dx;
  }

  return spacing;
}

}