// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef POLYTROPE_IC_H_PK11F
#define POLYTROPE_IC_H_PK11F

#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/local_equilibrium.hpp>

namespace zisa {

template <class EOS, class Gravity>
class PolytropeIC {
private:
  using eos_t = EOS;
  using gravity_t = Gravity;

public:
  PolytropeIC(std::shared_ptr<eos_t> eos, std::shared_ptr<gravity_t> gravity)
      : eos(std::move(eos)), gravity(std::move(gravity)) {}

  RhoP operator()(const XYZ &x) const {

    double alpha = gravity->alpha(zisa::norm(x));
    double eps = std::numeric_limits<double>::min();
    double r_eff = alpha * (zisa::norm(x) + eps);

    double rho = zisa::sin(r_eff) / r_eff;
    double p = zisa::pow<2>(rho);

    return RhoP{rho, p};
  }

private:
  std::shared_ptr<eos_t> eos;
  std::shared_ptr<gravity_t> gravity;
};

template <class EOS, class Gravity>
class GeneralPolytropeIC {
private:
  using eos_t = EOS;
  using gravity_t = Gravity;

public:
  GeneralPolytropeIC(const std::shared_ptr<EOS> &eos,
                     const std::shared_ptr<Gravity> &gravity,
                     const RhoEntropy &rhoK_center)
      : eos(eos), eq(eos, gravity), x_ref(XYZ::zeros()) {

    theta_ref = eos->enthalpy_entropy(rhoK_center);
  }

  RhoP operator()(const XYZ &x) const {
    return eos->rhoP(eq.extrapolate(theta_ref, x_ref, x));
  }

private:
  std::shared_ptr<eos_t> eos;
  IsentropicEquilibrium<eos_t, gravity_t> eq;
  EnthalpyEntropy theta_ref;
  XYZ x_ref;
};

template <class EOS, class Gravity>
class PolytropeWithJumpIC {
private:
  using eos_t = EOS;
  using gravity_t = Gravity;
  using eq_t = IsentropicEquilibrium<eos_t, gravity_t>;

public:
  PolytropeWithJumpIC(const std::shared_ptr<eos_t> &eos,
                      const std::shared_ptr<gravity_t> &gravity,
                      const EnthalpyEntropy &theta_inner,
                      const EnthalpyEntropy &theta_outer,
                      const XYZ &x_ref,
                      double amp_shape,
                      int n_bumps)
      : eos(eos),
        x_ref(x_ref),
        inner_equilibrium(eq_t(eos, gravity), theta_inner, x_ref),
        outer_equilibrium(eq_t(eos, gravity), theta_outer, x_ref),
        amp_shape(amp_shape),
        n_bumps(n_bumps) {}

  RhoP operator()(const XYZ &x) const {
    const auto &eq = (is_inner(x) ? inner_equilibrium : outer_equilibrium);
    RhoE rhoE = eq.extrapolate(x);
    return eos->rhoP(rhoE);
  }

  bool is_inner(const XYZ &x) const {
    double alpha = n_bumps * zisa::angle(x);
    double r_ref = zisa::norm(x_ref) + amp_shape * zisa::sin(alpha);

    return zisa::norm(x) < r_ref;
  }

private:
  std::shared_ptr<eos_t> eos;
  XYZ x_ref;
  LocalEquilibrium<eq_t> inner_equilibrium;
  LocalEquilibrium<eq_t> outer_equilibrium;

  double amp_shape;
  int n_bumps;
};

} // namespace zisa
#endif /* end of include guard */
