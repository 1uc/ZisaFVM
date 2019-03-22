#ifndef POLYTROPE_IC_H_PK11F
#define POLYTROPE_IC_H_PK11F

#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/local_equilibrium.hpp>

namespace zisa {

template <class EULER>
class PolytropeIC {
private:
  using euler_t = EULER;

public:
  PolytropeIC(const euler_t &euler) : euler(euler) {}

  RhoP operator()(const XYZ &x) const {

    double alpha = this->euler.gravity.alpha(zisa::norm(x));
    double eps = std::numeric_limits<double>::min();
    double r_eff = alpha * (zisa::norm(x) + eps);

    double rho = zisa::sin(r_eff) / r_eff;
    double p = zisa::pow<2>(rho);

    return RhoP{rho, p};
  }

private:
  euler_t euler;
};

template <class EULER>
class PolytropeWithJumpIC {
private:
  using euler_t = EULER;
  using eos_t = typename euler_t::eos_t;
  using gravity_t = typename euler_t::gravity_t;
  using eq_t = IsentropicEquilibrium<eos_t, gravity_t>;

public:
  PolytropeWithJumpIC(const euler_t &euler,
                      const EnthalpyEntropy &theta_inner,
                      const EnthalpyEntropy &theta_outer,
                      const XYZ &x_ref)
      : euler(euler),
        x_ref(x_ref),
        inner_equilibrium(eq_t(euler, 0), theta_inner, x_ref),
        outer_equilibrium(eq_t(euler, 0), theta_outer, x_ref) {}

  RhoP operator()(const XYZ &x) const {
    const auto &equilibrium
        = (zisa::norm(x) < zisa::norm(x_ref) ? inner_equilibrium
                                             : outer_equilibrium);
    RhoE rhoE = equilibrium.extrapolate(x);
    return euler.eos.rhoP(rhoE);
  }

private:
  euler_t euler;
  XYZ x_ref;
  LocalEquilibrium<eq_t> inner_equilibrium;
  LocalEquilibrium<eq_t> outer_equilibrium;
};

} // namespace zisa
#endif /* end of include guard */
