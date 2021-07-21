// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

/*
 *
 */

#include <zisa/math/poly2d.hpp>

namespace zisa {

ANY_DEVICE int_t poly_dof(int deg, int n_dims) {
  if (n_dims == 2) {
    return poly_dof<2>(deg);
  } else {
    return poly_dof<3>(deg);
  }
}

int poly_degree(int_t n_coeffs, int n_dims) {
  if (n_dims == 2) {
    return poly_degree<2>(n_coeffs);
  } else {
    return poly_degree<3>(n_coeffs);
  }
}

} // namespace zisa
