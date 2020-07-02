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
#include <zisa/model/equation_of_state.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

class HelmholtzEOS : public EquationOfState {
private:
  using super = EquationOfState;

public:
  using cvars_t = super::cvars_t;
  using xvars_t = super::xvars_t;

public:
  HelmholtzEOS() {}

  HelmholtzEOS(const std::string &) { LOG_ERR("Switch to other constructor."); }

  HelmholtzEOS(const std::string &eos_table,
               const std::vector<double> &mass_mixing_ratio,
               const std::vector<double> &mass_number,
               const std::vector<double> &charge_number) {
    int n_chars = int(eos_table.length());

    int status = -1;
    helmholtz_eos_readtable_c(eos_table.c_str(), &n_chars, &status);
    LOG_ERR_IF(
        status != 0,
        string_format("Reading EOS table failed. [%s]", eos_table.c_str()));

    int n_species = integer_cast<int>(charge_number.size());
    composition_bar_c(mass_mixing_ratio.data(),
                      mass_number.data(),
                      charge_number.data(),
                      &n_species,
                      &a_bar,
                      &z_bar,
                      &status);
    LOG_ERR_IF(status != 0,
               string_format("EOS failed to compute mean composition. [%s]",
                             eos_table.c_str()));
  }

  RhoE rhoE(const euler_var_t &u) const { return {u[0], internal_energy(u)}; }

  euler_var_t cvars(const RhoE &rhoE) const {
    return {rhoE.rho(), 0.0, 0.0, 0.0, rhoE.E()};
  }

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

    LOG_ERR_IF(tv1 <= 0, string_format("Negative `tv1`. [%e]", tv1));
    LOG_ERR_IF(tv2 <= 0, string_format("Negative `tv2`. [%e]", tv2));

    helmholtz_eos_bar_c(&id, &tv1, &tv2, &a_bar, &z_bar, &eos_ret, &status);
    LOG_ERR_IF(
        status != 0,
        string_format("EOS failed. %d [%d, %e, %e]", status, id, tv1, tv2));

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

    helmholtz_eos_wguess_bar_c(
        &TD_HS, &h, &s, &a_bar, &z_bar, &ret, &status, &rho_guess, &T_guess);

    LOG_ERR_IF(status != 0,
               string_format("EOS TD_HS failed. %d, [%.3e, %.3e; (%.3e, %.3e)]",
                             status,
                             h,
                             s,
                             rho_guess,
                             T_guess));

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
