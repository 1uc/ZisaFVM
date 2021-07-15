// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

namespace zisa {

WENOPoly WENO_AO::reconstruct(const array_view<double, 2, row_major> &rhs,
                              const array_view<WENOPoly, 1> &polys,
                              const array_const_view<cvars_t, 1> &qbar) const {

  auto shape = shape_t<2>(qbar.shape(0), cvars_t::size());
  auto ptr = (double *)(qbar.raw());
  auto view = array_const_view<double, 2>(shape, ptr);

  compute_polys(rhs, polys, view);
  return hybridize(array_const_view(polys));
}

ScalarPoly WENO_AO::reconstruct(const array_view<double, 2, row_major> &rhs,
                                const array_view<ScalarPoly, 1> &polys,
                                const array_const_view<double, 1> &qbar) const {

  auto shape = shape_t<2>(qbar.shape(0), 1);
  auto ptr = (double *)(qbar.raw());
  auto view = array_const_view<double, 2>(shape, ptr);

  compute_polys(rhs, polys, view);
  return hybridize(array_const_view(polys));
}

} // namespace zisa
