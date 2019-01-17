#ifndef LOCAL_EQUILIBRIUM_IMPL_H_NLD8L
#define LOCAL_EQUILIBRIUM_IMPL_H_NLD8L

#include "local_equilibrium_decl.hpp"
#include <zisa/math/quadrature.hpp>
#include <zisa/math/quasi_newton.hpp>

namespace zisa {

template <class Equilibrium>
RhoE extrapolate(const Equilibrium &eq,
                 const EnthalpyEntropy &theta,
                 const XY &xy_ref,
                 const XY &xy) {
  return extrapolate(eq, eq.eos, theta, xy_ref, xy);
}

template <class Equilibrium>
LocalEquilibrium<Equilibrium>::LocalEquilibrium(const Equilibrium &equilibrium,
                                                const Triangle &tri_ref)
    : tri_ref(tri_ref), equilibrium(equilibrium) {}

template <class Equilibrium>
void LocalEquilibrium<Equilibrium>::solve(const RhoE &rhoE_bar) {
  auto f = [this, &rhoE_bar](const EnthalpyEntropy &theta_star) {
    auto deg = equilibrium.quad_deg;

    auto rhoE_eq = [this, &theta_star](const XY &xy) {
      return zisa::extrapolate(
          equilibrium, theta_star, barycenter(tri_ref), xy);
    };

    auto vol = volume(tri_ref);
    auto rhoE_eq_bar
        = average(rhoE_eq, tri_ref, deg);

    return RhoE(rhoE_bar - rhoE_eq_bar);
  };

  auto df = [&f](const EnthalpyEntropy &x, int dir) {
    double eps = 1e-6 * zisa::abs(x[dir]);
    auto x_eps = static_cast<EnthalpyEntropy>(x + eps * x.unit_vector(dir));

    return RhoE((f(x_eps) - f(x)) / eps);
  };

  auto inv_df = [&df](const EnthalpyEntropy &x) {
    auto df0 = df(x, 0);
    auto df1 = df(x, 1);

    return [df0, df1](const auto &fx) {
      double inv_det = 1.0 / (df0[0] * df1[1] - df0[1] * df1[0]);
      return EnthalpyEntropy{inv_det * (df1[1] * fx[0] - df1[0] * fx[1]),
                             inv_det * (-df0[1] * fx[0] + df0[0] * fx[1])};
    };
  };

  const auto &eos = equilibrium.eos;
  auto guess = eos.enthalpy_entropy(rhoE_bar);

  auto atol = EnthalpyEntropy(1e-10 * guess);
  std::tie(theta, found_equilibrium) = quasi_newton(f, inv_df, guess, atol);
};

template <class Equilibrium>
RhoE LocalEquilibrium<Equilibrium>::extrapolate(const XY &xy) const {
  return zisa::extrapolate(equilibrium, theta, barycenter(tri_ref), xy);
}

template <class Equilibrium>
RhoE LocalEquilibrium<Equilibrium>::extrapolate(const Triangle &tri) const {
  auto deg = equilibrium.quad_deg;
  return average([this](const XY &xy) { return extrapolate(xy); }, tri, deg);
}
} // namespace zisa

#endif /* end of include guard */
