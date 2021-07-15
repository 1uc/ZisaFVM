// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>

namespace zisa {

WENOPoly CWENO_AO::reconstruct(const array_view<double, 2, row_major> &rhs,
                               const array_view<WENOPoly, 1> &polys,
                               const array_const_view<cvars_t, 1> &qbar) const {
  auto shape = shape_t<2>(qbar.shape(0), cvars_t::size());
  auto ptr = (double *)(qbar.raw());
  auto view = array_const_view<double, 2>(shape, ptr);

  return reconstruct_impl(rhs, polys, view);
}

ScalarPoly
CWENO_AO::reconstruct(const array_view<double, 2, row_major> &rhs,
                      const array_view<ScalarPoly, 1> &polys,
                      const array_const_view<double, 1> &qbar) const {
  auto shape = shape_t<2>(qbar.shape(0), 1);
  auto ptr = (double *)(qbar.raw());
  auto view = array_const_view<double, 2>(shape, ptr);

  return reconstruct_impl(rhs, polys, view);
}

template <class Poly>
Poly CWENO_AO::reconstruct_impl(const array_view<double, 2, row_major> &rhs,
                                const array_view<Poly, 1> &polys,
                                const array_const_view<double, 2> &qbar) const {
  compute_polys(rhs, polys, qbar);

  auto k_high = stencils.highest_order_stencil();
  for (int_t k = 0; k < stencils.size(); ++k) {
    if (k_high != k) {
      polys[k_high] -= linear_weights[k] * polys[k];
    }
  }

  // p_opt = C * p_c + \sum_i c_i p_i
  // p_c = (p_opt - \sum_i c_i p_i) / C
  polys[k_high] /= linear_weights[k_high];

  return hybridize(polys);
}

} // namespace zisa
