/*
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2018-01-24
 */

#ifndef HELMHOLTZ_EOS_H_BJMBP
#define HELMHOLTZ_EOS_H_BJMBP

#if ZISA_HAS_HELMHOLTZ_EOS != 1
#error "Missing `ZISA_HAS_HELMHOLTZ_EOS=1`."
#endif

#include "helmholtzeos_module_c.h"
#include <zisa/memory/array_view.hpp>
#include <zisa/model/equation_of_state.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

void initialize_helmholtz_eos(std::string &eos_table);

class HelmholtzEOS : public EquationOfState {
private:
  using super = EquationOfState;

public:
  using cvars_t = super::cvars_t;
  using xvars_t = super::xvars_t;

public:
  HelmholtzEOS() : a_bar(-1.0), z_bar(-1.0) {}

  HelmholtzEOS(double a_bar, double z_bar) : a_bar(a_bar), z_bar(z_bar) {}

  HelmholtzEOS(const array_const_view<double, 1> &mass_mixing_ratio,
               const array_const_view<double, 1> &mass_number,
               const array_const_view<double, 1> &charge_number) {

    int status = -1;
    int n_species = integer_cast<int>(charge_number.size());
    composition_bar_c(raw_ptr(mass_mixing_ratio),
                      raw_ptr(mass_number),
                      raw_ptr(charge_number),
                      &n_species,
                      &a_bar,
                      &z_bar,
                      &status);

    LOG_ERR_IF(status != 0 || !ispositive(a_bar) || !ispositive(z_bar),
               string_format("EOS failed: status = %d, a_bar = %e, z_bar = %e, "
                             "X = %s, A = %s, Z = %s",
                             status,
                             a_bar,
                             z_bar,
                             format_as_list(mass_mixing_ratio).c_str(),
                             format_as_list(mass_number).c_str(),
                             format_as_list(charge_number).c_str()));
  }

  using super::rhoE;
  RhoE rhoE(const euler_var_t &u) const { return {u[0], internal_energy(u)}; }

  euler_var_t cvars(const RhoE &rhoE) const {
    return {rhoE.rho(), 0.0, 0.0, 0.0, rhoE.E()};
  }

  using super::xvars;
  xvars_t xvars(const RhoE &rhoE) const { return xvars(cvars(rhoE)); }

  xvars_t xvars(const cvars_t &u) const {
    auto tmp = full_extra_variables(u);
    return xvars_t{tmp.p, tmp.a};
  }

  double pressure(const RhoE &rhoE) const {
    auto tmp = full_extra_variables(rhoE);
    return tmp.p;
  }

  RhoE rhoE(const EnthalpyEntropy &theta, const RhoT &rhoT_guess) const {
    auto xvars = full_extra_variables(theta, rhoT_guess);
    return RhoE{xvars.rho, xvars.E};
  }

  euler_full_xvars_t full_extra_variables(const euler_var_t &u) const {
    auto E_int = internal_energy(u);
    return full_extra_variables(RhoE{u[0], E_int});
  }

  euler_full_xvars_t full_extra_variables(const RhoE &rhoE) const {
    return full_extra_variables(TD_ED, rhoE.E() / rhoE.rho(), rhoE.rho());
  }

  euler_full_xvars_t full_extra_variables(const RhoP &rhoP) const {
    return full_extra_variables(TD_PD, rhoP.p(), rhoP.rho());
  }

  euler_full_xvars_t full_extra_variables(const RhoEntropy &rhoS) const {
    return full_extra_variables(TD_SD, rhoS.s(), rhoS.rho());
  }

  euler_full_xvars_t full_extra_variables(const RhoT &rhoT) const {
    return full_extra_variables(TD_TD, rhoT.T(), rhoT.rho());
  }

  euler_full_xvars_t
  full_extra_variables(int id, double tv1, double tv2) const {
    int status = -1;
    eos_state_type eos_ret;

    LOG_ERR_IF(
        (tv1 <= 0.0) || (tv2 <= 0.0),
        string_format(
            "Invalid input to EOS. %d : tv1 = %e, tv2 = %e", id, tv1, tv2));

    helmholtz_eos_bar_c(&id, &tv1, &tv2, &a_bar, &z_bar, &eos_ret, &status);
    LOG_ERR_IF(
        status != 0,
        string_format("EOS failed. %d [%d, %e, %e] a_bar = %e, z_bar = %e",
                      status,
                      id,
                      tv1,
                      tv2,
                      a_bar,
                      z_bar));

    return full_extra_variables(eos_ret);
  }

  euler_full_xvars_t full_extra_variables(const EnthalpyEntropy &theta,
                                          const RhoT &rhoT_guess) const {

    int status = -1;
    eos_state_type ret;

    double h = theta.h();
    double s = theta.s();

    double rho_guess = rhoT_guess.rho();
    double T_guess = rhoT_guess.T();

    LOG_ERR_IF(
        ((h <= 0.0) || (s <= 0.0) || (rho_guess <= 0.0) || (T_guess <= 0.0)),
        string_format("Invalid input to EOS. h = %e, s = %e; rho_guess "
                      "= %e, T_guess = %e",
                      h,
                      s,
                      rho_guess,
                      T_guess));

    helmholtz_eos_wguess_bar_c(
        &TD_HS, &h, &s, &a_bar, &z_bar, &ret, &status, &T_guess, &rho_guess);

    LOG_ERR_IF(status != 0,
               string_format("EOS TD_HS failed. %d, [%.3e, %.3e; (%.3e, %.3e)] "
                             "a_bar = %e, z_bar = %e",
                             status,
                             h,
                             s,
                             rho_guess,
                             T_guess,
                             a_bar,
                             z_bar));

    return full_extra_variables(ret);
  }

  euler_full_xvars_t full_extra_variables(const eos_state_type &eos_ret) const {

    euler_full_xvars_t ret;
    ret.rho = eos_ret.d;
    ret.E = eos_ret.d * eos_ret.e;
    ret.a = eos_ret.cs;
    ret.h = eos_ret.h;
    ret.s = eos_ret.s;
    ret.p = eos_ret.p;
    ret.T = eos_ret.t;

    return ret;
  }

  std::string str() const { return "Helmholtz EOS."; }

private:
  double a_bar;
  double z_bar;
};

void save(HDF5Writer &writer, const HelmholtzEOS &eos);

} // namespace zisa
#endif /* end of include guard */
