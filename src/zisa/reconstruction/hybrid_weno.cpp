#include <algorithm>
#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/reconstruction/hybrid_weno.hpp>

namespace zisa {

HybridWENO::HybridWENO(const std::shared_ptr<Grid> &grid,
                       int_t i_cell,
                       const HybridWENO_Params &params)
    : stencils(grid, i_cell, params.stencil_family_params),
      lsq_solvers(grid, stencils),
      polys(stencils.size()),
      linear_weights(stencils.size()),
      epsilon(1e-6),
      exponent(4) {
  rhs = array<double, 1>(shape_t<1>{stencils.combined_stencil_size()});

  auto tot = std::accumulate(
      params.linear_weights.begin(), params.linear_weights.end(), 0.0);
  for (int_t i = 0; i < linear_weights.size(); ++i) {
    linear_weights[i] = params.linear_weights[i] / tot;
  }
}

void HybridWENO::compute_polys(const array<double, 1> &qbar) const {

  auto p_avg = Poly2D<MAX_DEGREE>{{qbar(0)}, {0.0}};

  for (int_t k = 0; k < stencils.size(); ++k) {
    for (int_t i = 0; i < stencils[k].size() - 1; ++i) {
      rhs(i) = qbar(stencils[k].local(i + 1)) - qbar(0);
    }

    polys[k] = p_avg + lsq_solvers[k].solve(rhs);
  }
}

auto HybridWENO::hybridize() const -> Poly2D<MAX_DEGREE> {
  return eno_hybridize();
}

auto HybridWENO::eno_hybridize() const -> Poly2D<MAX_DEGREE> {
  double al_tot = 0.0;
  auto p = Poly2D<MAX_DEGREE>{{0.0}, {0.0}};
  for (int_t k = 0; k < stencils.size(); ++k) {
    auto IS = smoothness_indicator(polys[k]);

    auto al = linear_weights[k] / (epsilon + zisa::pow(IS, exponent));
    al_tot += al;

    p += al * polys[k];
  }

  return p / al_tot;
}

auto HybridWENO::tau_hybridize() const -> Poly2D<MAX_DEGREE> {
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

auto HybridWENO::local2global() const -> decltype(stencils.local2global()) {
  return stencils.local2global();
}

} // namespace zisa
