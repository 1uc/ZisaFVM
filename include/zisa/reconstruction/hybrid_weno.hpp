// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

/* Base class for hybridized WENO.
 */

#ifndef HYBRID_WENO_H_49M5Q
#define HYBRID_WENO_H_49M5Q

#include <zisa/config.hpp>

#include <zisa/math/poly2d.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/lsq_solver_family.hpp>
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {
class HybridWENO {
protected:
  using cvars_t = euler_var_t;

protected:
  StencilFamily stencils;
  LSQSolverFamily lsq_solvers;

  array<double, 1> linear_weights;
  mutable array<double, 1> non_linear_weights;
  double epsilon;
  double exponent;

public:
  HybridWENO() = default;

  HybridWENO(const std::shared_ptr<Grid> &grid,
             int_t i_cell,
             const HybridWENOParams &params);

  HybridWENO(const std::shared_ptr<Grid> &grid,
             StencilFamily stencil_family,
             int_t i_cell,
             const HybridWENOParams &params);

  const std::vector<int_t> &local2global() const;
  int_t combined_stencil_size() const;

  /// Indistinguishable by calls to the public interface.
  bool operator==(const HybridWENO &other) const;

  /// Can be distinguished by calls to the public interface.
  bool operator!=(const HybridWENO &other) const;

  std::string str(int /* verbose */ = 0) const {
    return string_format("%s, eps = %e, s = %f",
                         format_as_list(linear_weights).c_str(),
                         epsilon,
                         exponent);
  }

protected:
  WENOPoly hybridize(const array_const_view<WENOPoly, 1> &polys) const;
  ScalarPoly hybridize(const array_const_view<ScalarPoly, 1> &polys) const;

  void compute_polys(const array_view<double, 2, row_major> &rhs,
                     const array_view<WENOPoly, 1> &polys,
                     const array_const_view<double, 2> &qbar) const;

  void compute_polys(const array_view<double, 2, row_major> &rhs,
                     const array_view<ScalarPoly, 1> &polys,
                     const array_const_view<double, 2> &qbar) const;

private:
  template <class Poly>
  Poly hybridize_impl(const array_const_view<Poly, 1> &polys) const;

  template <class Poly>
  void compute_polys_impl(const array_view<double, 2, row_major> &rhs,
                          const array_view<Poly, 1> &polys,
                          const array_const_view<double, 2> &qbar) const;

  template <class Poly>
  Poly eno_hybridize(const array_const_view<Poly, 1> &polys) const;

  template <class Poly>
  Poly tau_hybridize(array<Poly, 1> &polys) const;
};

} // namespace zisa
#endif /* end of include guard */
