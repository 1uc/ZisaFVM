#ifndef EQUATION_OF_STATE_BASE_H_F1NYM
#define EQUATION_OF_STATE_BASE_H_F1NYM

#include <zisa/config.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

class EquationOfState {
public:
  using cvars_t = euler_var_t;
  using xvars_t = cvars_t::xvars_t;

public:
  RhoE rhoE(const euler_full_xvars_t &full_xvars) const {
    return {full_xvars.rho, full_xvars.E};
  }

  xvars_t xvars(const euler_full_xvars_t &full_xvars) const {
    auto xvars_ = xvars_t{};
    xvars_.p = full_xvars.p;
    xvars_.a = full_xvars.a;
    return xvars_;
  }

  double kinetic_energy(const cvars_t &u) const {
    return zisa::kinetic_energy(u);
  }

  double internal_energy(const cvars_t &u) const {
    return zisa::internal_energy(u);
  }
};

} // namespace zisa

#endif /* end of include guard */
