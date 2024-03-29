// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <algorithm>
#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/reconstruction/hybrid_weno.hpp>

namespace zisa {

HybridWENO::HybridWENO(const std::shared_ptr<Grid> &grid,
                       int_t i_cell,
                       const HybridWENOParams &params)
    : stencils(*grid, i_cell, params.stencil_family_params),
      lsq_solvers(grid, stencils),
      linear_weights(stencils.size()),
      non_linear_weights(stencils.size()),
      epsilon(params.epsilon),
      exponent(params.exponent) {
  assert(params.linear_weights.size() == stencils.size());

  auto tot = std::accumulate(
      params.linear_weights.begin(), params.linear_weights.end(), 0.0);

  for (int_t i = 0; i < linear_weights.size(); ++i) {
    linear_weights[i] = params.linear_weights[i] / tot;
  }
}

HybridWENO::HybridWENO(const std::shared_ptr<Grid> &grid,
                       StencilFamily stencil_family,
                       int_t /* i_cell */,
                       const HybridWENOParams &params)
    : stencils(std::move(stencil_family)),
      lsq_solvers(grid, stencils),
      linear_weights(stencils.size()),
      non_linear_weights(stencils.size()),
      epsilon(params.epsilon),
      exponent(params.exponent) {
  assert(params.linear_weights.size() == stencils.size());

  auto tot = std::accumulate(
      params.linear_weights.begin(), params.linear_weights.end(), 0.0);

  for (int_t i = 0; i < linear_weights.size(); ++i) {
    linear_weights[i] = params.linear_weights[i] / tot;
  }
}

int_t HybridWENO::combined_stencil_size() const {
  return stencils.combined_stencil_size();
}

void HybridWENO::compute_polys(const array_view<double, 2, row_major> &rhs,
                               const array_view<WENOPoly, 1> &polys,
                               const array_const_view<double, 2> &qbar) const {
  return compute_polys_impl<WENOPoly>(rhs, polys, qbar);
}

void HybridWENO::compute_polys(const array_view<double, 2, row_major> &rhs,
                               const array_view<ScalarPoly, 1> &polys,
                               const array_const_view<double, 2> &qbar) const {

  return compute_polys_impl<ScalarPoly>(rhs, polys, qbar);
}

template <class Poly>
void HybridWENO::compute_polys_impl(
    const array_view<double, 2, row_major> &rhs,
    const array_view<Poly, 1> &polys,
    const array_const_view<double, 2> &qbar) const {

  for (int_t k = 0; k < stencils.size(); ++k) {
    for (int_t ig = 0; ig < stencils[k].size() - 1; ++ig) {
      int_t il = stencils[k].local(ig + 1);

      for (int_t k_var = 0; k_var < Poly::n_vars(); ++k_var) {
        rhs(ig, k_var) = qbar(il, k_var) - qbar(0, k_var);
      }
    }

    polys[k] = lsq_solvers[k].solve<Poly>(rhs);

    for (int_t k_var = 0; k_var < Poly::n_vars(); ++k_var) {
      polys[k].a(k_var) = qbar(0, k_var);
    }
  }
}

WENOPoly
HybridWENO::hybridize(const array_const_view<WENOPoly, 1> &polys) const {
  return hybridize_impl<WENOPoly>(polys);
}

ScalarPoly
HybridWENO::hybridize(const array_const_view<ScalarPoly, 1> &polys) const {
  return hybridize_impl<ScalarPoly>(polys);
}

template <class Poly>
Poly HybridWENO::hybridize_impl(const array_const_view<Poly, 1> &polys) const {
  return eno_hybridize<Poly>(polys);
}

template <class Poly>
Poly HybridWENO::eno_hybridize(const array_const_view<Poly, 1> &polys) const {
  double al_tot = 0.0;
  for (int_t k = 0; k < stencils.size(); ++k) {
    auto IS = zisa::maximum(smoothness_indicator(polys[k]));

    double al = linear_weights[k] / (epsilon + zisa::pow(IS, exponent));

    non_linear_weights[k] = al;
    al_tot += al;
  }

  int n_dims = 2;
  auto p = Poly(0, {0.0}, XYZ::zeros(), 1.0, n_dims);
  for (int_t k = 0; k < stencils.size(); ++k) {
    p += (non_linear_weights[k] / al_tot) * polys[k];
  }

  return p;
}

template <class Poly>
Poly HybridWENO::tau_hybridize(array<Poly, 1> &) const {
  LOG_ERR("Not implemented.");
  // auto k_high = stencils.highest_order_stencil();

  // auto beta_high = smoothness_indicator(polys[k_high]);
  // for (int_t k = 0; k < stencils.size(); ++k) {
  //   if (k != k_high) {
  //     auto beta = smoothness_indicator(polys[k]);
  //     tau += zisa::abs(beta_high - beta);
  //   }
  // }

  // tau /= stencils.size() - 1;
}

bool HybridWENO::operator==(const HybridWENO &other) const {
  if (epsilon != other.epsilon) {
    return false;
  }
  if (exponent != other.exponent) {
    return false;
  }
  if (linear_weights != other.linear_weights) {
    return false;
  }
  if (stencils != other.stencils) {
    return false;
  }

  return lsq_solvers == other.lsq_solvers;
}

bool HybridWENO::operator!=(const HybridWENO &other) const {
  return !((*this) == other);
}

const std::vector<int_t> &HybridWENO::local2global() const {
  return stencils.local2global();
}

} // namespace zisa
