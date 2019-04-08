#ifndef LOCAL_EQUILIBRIUM_IMPL_H_NLD8L
#define LOCAL_EQUILIBRIUM_IMPL_H_NLD8L

#include "local_equilibrium_decl.hpp"
#include <zisa/math/quadrature.hpp>
#include <zisa/math/quasi_newton.hpp>

namespace zisa {

template <class Equilibrium>
LocalEquilibrium<Equilibrium>::LocalEquilibrium(const Equilibrium &equilibrium)
    : equilibrium(equilibrium) {}

template <class Equilibrium>
LocalEquilibrium<Equilibrium>::LocalEquilibrium(
    const Equilibrium &equilibrium,
    const EnthalpyEntropy &theta_ref,
    const XYZ &x_ref)
    : theta(theta_ref),
      x_ref(x_ref),
      found_equilibrium(true),
      equilibrium(equilibrium) {}

template <class Equilibrium>
void LocalEquilibrium<Equilibrium>::solve(const RhoE &rhoE_bar,
                                          const Triangle &tri_ref) {

  x_ref = barycenter(tri_ref);

  auto f = [this, &tri_ref, &rhoE_bar](const EnthalpyEntropy &theta_star) {
    auto deg = equilibrium.quad_deg;

    auto rhoE_eq = [this, &theta_star](const XYZ &xy) {
      return equilibrium.extrapolate(theta_star, x_ref, xy);
    };

    return RhoE(rhoE_bar - average(rhoE_eq, tri_ref, deg));
  };

  auto df = [&f](const EnthalpyEntropy &x, int_t dir) {
    double eps = 1e-6 * zisa::abs(x[dir]);
    auto x_eps = static_cast<EnthalpyEntropy>(x + eps * x.unit_vector(dir));

    return RhoE((f(x_eps) - f(x)) / eps);
  };

  auto inv_df = [&df](const EnthalpyEntropy &x) {
    auto df0 = df(x, 0ul);
    auto df1 = df(x, 1ul);

    return [df0, df1](const auto &fx) {
      double inv_det = 1.0 / (df0[0] * df1[1] - df0[1] * df1[0]);
      return EnthalpyEntropy{inv_det * (df1[1] * fx[0] - df1[0] * fx[1]),
                             inv_det * (-df0[1] * fx[0] + df0[0] * fx[1])};
    };
  };

  const auto &eos = equilibrium.euler->eos;
  auto guess = eos.enthalpy_entropy(rhoE_bar);

  auto atol = EnthalpyEntropy(1e-10 * guess);
  std::tie(theta, found_equilibrium) = quasi_newton(f, inv_df, guess, atol);
};

template <class Equilibrium>
RhoE LocalEquilibrium<Equilibrium>::extrapolate(const XYZ &xy) const {
  return equilibrium.extrapolate(theta, x_ref, xy);
}

template <class Equilibrium>
RhoE LocalEquilibrium<Equilibrium>::extrapolate(const Triangle &tri) const {
  auto deg = equilibrium.quad_deg;
  return average([this](const XYZ &xy) { return extrapolate(xy); }, tri, deg);
}

} // namespace zisa

#endif /* end of include guard */
