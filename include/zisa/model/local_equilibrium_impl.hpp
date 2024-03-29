// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef LOCAL_EQUILIBRIUM_IMPL_H_NLD8L
#define LOCAL_EQUILIBRIUM_IMPL_H_NLD8L

#include "local_equilibrium_decl.hpp"

#include <zisa/math/cell.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/math/quasi_newton.hpp>

#include <functional>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/src/NumericalDiff/NumericalDiff.h>

namespace zisa {

template <class Equilibrium>
LocalEquilibriumBase<Equilibrium>::LocalEquilibriumBase(
    const Equilibrium &equilibrium)
    : equilibrium(equilibrium) {}

template <class Equilibrium>
LocalEquilibriumBase<Equilibrium>::LocalEquilibriumBase(
    const Equilibrium &equilibrium,
    const equilibrium_values_t &theta_ref,
    const XYZ &x_ref)
    : theta(theta_ref),
      x_ref(x_ref),
      found_equilibrium(true),
      equilibrium(equilibrium) {}

template <class Equilibrium>
void LocalEquilibriumBase<Equilibrium>::solve(const RhoE &rhoE_bar,
                                              const Cell &cell_ref) {
  solve_exact(rhoE_bar, cell_ref);
}

template <class Equilibrium>
void LocalEquilibriumBase<Equilibrium>::solve_exact(const RhoE &rhoE_bar,
                                                    const Cell &cell_ref) {

  x_ref = cell_ref.qr.points[0];

  const auto &eos = *equilibrium.eos;

  auto full_guess = eos.full_extra_variables(rhoE_bar);
  auto enthalpy_entropy_guess = EnthalpyEntropy{full_guess.h, full_guess.s};
  auto rhoT_guess = RhoT{full_guess.rho, full_guess.T};

  auto f = [this, &eos, &cell_ref, &rhoE_bar, &rhoT_guess](
               const EnthalpyEntropy &theta_star) {
    auto rhoE_eq = [this, &eos, &theta_star, &rhoT_guess](const XYZ &xy) {
      return equilibrium.extrapolate(
          isentropic_equilibrium_values(eos, theta_star, rhoT_guess),
          x_ref,
          xy);
    };

    return RhoE(rhoE_bar - average(cell_ref, rhoE_eq));
  };

  auto df = [&f](const EnthalpyEntropy &x, int_t dir) {
    double eps = 1e-6 * zisa::abs(x[dir]);

    auto x_p_eps
        = static_cast<EnthalpyEntropy>(x + 0.5 * eps * x.unit_vector(dir));

    auto x_m_eps
        = static_cast<EnthalpyEntropy>(x - 0.5 * eps * x.unit_vector(dir));

    return RhoE((f(x_p_eps) - f(x_m_eps)) / eps);
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

  auto atol = EnthalpyEntropy(1e-13 * enthalpy_entropy_guess);

  auto [hS, has_eq] = quasi_newton(f, inv_df, enthalpy_entropy_guess, atol);

  found_equilibrium = has_eq;
  theta = isentropic_equilibrium_values(eos, hS, rhoT_guess);
}

namespace detail {

class LMDRhoEBar : public Eigen::DenseFunctor<double> {
public:
  LMDRhoEBar(const std::function<RhoE(const EnthalpyEntropy &)> &f_)
      : DenseFunctor<double>(2, 2), f(f_) {}

  int operator()(const Eigen::VectorXd &theta,
                 Eigen::VectorXd &drhoE_bar) const {

    auto drhoE_bar_ = f(EnthalpyEntropy{theta[0], theta[1]});

    drhoE_bar[0] = drhoE_bar_.rho();
    drhoE_bar[1] = drhoE_bar_.E();

    return 0;
  }

private:
  std::function<RhoE(const EnthalpyEntropy &)> f;
};

}

template <class Equilibrium>
void LocalEquilibriumBase<Equilibrium>::solve_lsq(const RhoE &rhoE_bar,
                                                  const Cell &cell_ref) {

  x_ref = cell_ref.qr.points[0];

  const auto &eos = *equilibrium.eos;

  auto full_guess = eos.full_extra_variables(rhoE_bar);
  auto rhoT_guess = RhoT{full_guess.rho, full_guess.T};

  auto drhoE_bar = [this, &eos, &cell_ref, &rhoE_bar, &rhoT_guess](
                       const EnthalpyEntropy &theta_star) {
    auto rhoE_eq = [this, &eos, &theta_star, &rhoT_guess](const XYZ &xy) {
      return equilibrium.extrapolate(
          isentropic_equilibrium_values(eos, theta_star, rhoT_guess),
          x_ref,
          xy);
    };

    return RhoE(rhoE_bar - average(cell_ref, rhoE_eq));
  };

  auto lm_drhoE_bar = detail::LMDRhoEBar(drhoE_bar);
  auto lm_drhoE_bar_der
      = Eigen::NumericalDiff<detail::LMDRhoEBar>(lm_drhoE_bar);
  auto lm = Eigen::LevenbergMarquardt<Eigen::NumericalDiff<detail::LMDRhoEBar>>(
      lm_drhoE_bar_der);

  auto hS = Eigen::VectorXd(2);
  hS[0] = full_guess.h;
  hS[1] = full_guess.s;
  auto info = lm.minimize(hS);
  LOG_ERR_IF(info != 1, "`LevenbergMarquardt` failed.");

  found_equilibrium = true;
  theta = isentropic_equilibrium_values(
      eos, EnthalpyEntropy{hS[0], hS[1]}, rhoT_guess);
}

template <class Equilibrium>
RhoE LocalEquilibriumBase<Equilibrium>::extrapolate(const XYZ &xy) const {
  if (found_equilibrium) {
    return equilibrium.extrapolate(theta, x_ref, xy);
  } else {
    return RhoE(RhoE::zeros());
  }
}

template <class Equilibrium>
std::pair<RhoE, euler_xvar_t>
LocalEquilibriumBase<Equilibrium>::extrapolate_full(const XYZ &xy) const {
  if (found_equilibrium) {
    return equilibrium.extrapolate_full(theta, x_ref, xy);
  } else {
    return {RhoE(RhoE::zeros()), xvars_t{}};
  }
}

template <class Equilibrium>
RhoE LocalEquilibriumBase<Equilibrium>::extrapolate(const Cell &cell) const {
  if (found_equilibrium) {
    return average(cell, [this](const XYZ &xy) { return extrapolate(xy); });
  } else {
    return RhoE(RhoE::zeros());
  }
}

template <class Equilibrium>
std::string LocalEquilibriumBase<Equilibrium>::str(int verbose) const {
  if (verbose == 0) {
    if (found_equilibrium) {
      return "No equilibrium found.";
    } else {
      return string_format("theta = %s @ x_ref = %s",
                           format_as_list(theta).c_str(),
                           format_as_list(x_ref).c_str());
    }
  } else if (verbose >= 1) {
    return string_format("theta = %s @ x_ref = %s; found_equilibrium = %d",
                         format_as_list(theta).c_str(),
                         format_as_list(x_ref).c_str(),
                         found_equilibrium);
  }

  return "";
}

} // namespace zisa

#endif /* end of include guard */
