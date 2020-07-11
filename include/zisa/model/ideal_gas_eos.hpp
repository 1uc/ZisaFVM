/* Ideal Gas equation of state.
 */

#ifndef IDEAL_GAS_EOS_H_SNFCI
#define IDEAL_GAS_EOS_H_SNFCI

#include <zisa/model/equation_of_state_base.hpp>

namespace zisa {

class IdealGasEOS : public EquationOfState {
public:
  using cvars_t = euler_var_t;
  using xvars_t = cvars_t::xvars_t;

public:
  IdealGasEOS() = default;
  IdealGasEOS(double gamma, double specific_gas_constant)
      : gamma_(gamma), specific_gas_constant_(specific_gas_constant) {}

  xvars_t xvars(const RhoE &rhoE) const { return xvars(cvars(rhoE)); }

  xvars_t xvars(const cvars_t &u) const {
    xvars_t xvars;

    xvars.p = pressure(u);
    xvars.a = sound_speed(u);

    return xvars;
  }

  ANY_DEVICE_INLINE double gamma() const { return gamma_; }

  /// This is `kappa`.
  /** Sorry I don't think this has a proper name. */
  ANY_DEVICE_INLINE double kappa(void) const {
    return (gamma() - 1.0) / gamma();
  }

  ANY_DEVICE_INLINE double specific_gas_constant() const {
    return specific_gas_constant_;
  }

  /// Ratio of specific heat at constant pressure.
  ANY_DEVICE_INLINE double cP(void) const {
    return specific_gas_constant() / kappa();
  }

  /// Ratio of specific heat at constant volume.
  ANY_DEVICE_INLINE double cV(void) const {
    return 1.0 / (gamma() - 1.0) * specific_gas_constant();
  }

  ANY_DEVICE_INLINE euler_var_t cvars(RhoP rhoP) const {
    return cvars_t{rhoP.rho(), 0.0, 0.0, 0.0, internal_energy(rhoP)};
  }

  ANY_DEVICE_INLINE euler_var_t cvars(RhoE rhoE) const {
    const auto &[rho, E] = rhoE;
    return cvars_t{rho, 0.0, 0.0, 0.0, E};
  }

  ANY_DEVICE_INLINE EnthalpyEntropy enthalpy_entropy(RhoEntropy rhoK) const {
    return enthalpy_entropy(rhoE(rhoK));
  }

  ANY_DEVICE_INLINE EnthalpyEntropy enthalpy_entropy(RhoP rhoP) const {
    return EnthalpyEntropy{enthalpy(rhoP), entropy(rhoP)};
  }

  ANY_DEVICE_INLINE EnthalpyEntropy enthalpy_entropy(RhoE rhoE) const {
    return enthalpy_entropy(rhoP(rhoE));
  }

  ANY_DEVICE_INLINE RhoE rhoE(const euler_var_t &u) const {
    return {u[0], internal_energy(u)};
  }

  ANY_DEVICE_INLINE RhoE rhoE(EnthalpyEntropy theta, RhoT) const {
    return rhoE(theta);
  }

  ANY_DEVICE_INLINE RhoE rhoE(EnthalpyEntropy theta) const {
    auto rho = this->rho(theta);
    auto E = this->internal_energy(theta);

    return RhoE{rho, E};
  }

  ANY_DEVICE_INLINE RhoE rhoE(PressureEntropy theta, RhoT) const {
    return rhoE(theta);
  }

  ANY_DEVICE_INLINE RhoE rhoE(const PressureEntropy &theta) const {
    auto rho = this->rho(theta);
    auto E = this->internal_energy(theta);

    return RhoE{rho, E};
  }

  ANY_DEVICE_INLINE RhoE rhoE(const RhoP &rhoP) const {
    return RhoE{rhoP.rho(), internal_energy(rhoP)};
  }

  ANY_DEVICE_INLINE RhoE rhoE(const RhoEntropy &rhoK) const {
    double rho = rhoK.rho();
    return RhoE{rho, internal_energy(rhoK)};
  }

  ANY_DEVICE_INLINE RhoT rhoT(RhoE rhoE) const {
    auto rho = rhoE.rho();
    auto T = this->temperature(rhoE);

    return RhoT{rho, T};
  }

  ANY_DEVICE_INLINE RhoT rhoT(RhoEntropy rhoK) const {
    auto rho = rhoK.rho();
    auto T = this->temperature(rhoP(rhoK));

    return RhoT{rho, T};
  }

  ANY_DEVICE_INLINE RhoP rhoP(const euler_var_t &u) const {
    return RhoP{u[0], pressure(u)};
  }

  ANY_DEVICE_INLINE RhoP rhoP(RhoE rhoE) const {
    return RhoP{rhoE.rho(), pressure(rhoE)};
  }

  ANY_DEVICE_INLINE RhoP rhoP(RhoEntropy rhoK) const {
    return RhoP{rhoK.rho(), pressure(rhoK)};
  }

  ANY_DEVICE_INLINE double entropy(RhoP rhoP) const { return K(rhoP); }
  ANY_DEVICE_INLINE double entropy(RhoE rhoE) const { return K(rhoE); }

  ANY_DEVICE_INLINE double enthalpy(const RhoEntropy &rhoK) const {
    return enthalpy(rhoP(rhoK));
  }

  ANY_DEVICE_INLINE double enthalpy(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    return gamma() / (gamma() - 1.0) * p / rho;
  }

  ANY_DEVICE_INLINE double enthalpy(const RhoE &rhoE) const {
    return enthalpy(rhoP(rhoE));
  }

  ANY_DEVICE_INLINE double K(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    return p / zisa::pow(rho, gamma());
  }

  ANY_DEVICE_INLINE double K(const RhoE &rhoE) const { return K(rhoP(rhoE)); }

  ANY_DEVICE_INLINE double rho(EnthalpyEntropy theta) const {
    const auto &[h, K] = theta;
    double gamma = IdealGasEOS::gamma();
    double base = 1.0 / K * (gamma - 1.0) / gamma * h;
    double exponent = 1.0 / (gamma - 1.0);
    return zisa::pow(base, exponent);
  }

  ANY_DEVICE_INLINE double rho(PressureEntropy theta) const {
    const auto &[p, K] = theta;
    return zisa::pow(p / K, 1.0 / gamma());
  }

  /// Pressure given density and internal energy, using ideal gas law.

  ANY_DEVICE_INLINE double pressure(EnthalpyEntropy theta) const {
    double K = theta.s();
    double rho = this->rho(theta);

    return pressure(RhoEntropy{rho, K});
  }

  ANY_DEVICE_INLINE double pressure(RhoE rhoE) const {
    return rhoE.E() * (gamma() - 1.0);
  }

  ANY_DEVICE_INLINE double pressure(RhoEntropy rhoK) const {
    const auto &[rho, K] = rhoK;
    return K * zisa::pow(rho, gamma());
  }

  ANY_DEVICE_INLINE double pressure(RhoT rhoT) const {
    const auto &[rho, T] = rhoT;
    double r_gas = specific_gas_constant();
    return r_gas * rho * T;
  }

  ANY_DEVICE_INLINE double pressure(const euler_var_t &u) const {
    double E_kin = kinetic_energy(u);
    return pressure(RhoE{u[0], u[4] - E_kin});
  }

  using EquationOfState::internal_energy;
  using EquationOfState::kinetic_energy;

  /// Internal energy.
  ANY_DEVICE_INLINE double ideal_internal_energy(double p) const {
    return p / (gamma() - 1.0);
  }

  ANY_DEVICE_INLINE double internal_energy(EnthalpyEntropy theta) const {
    double p = pressure(theta);
    return ideal_internal_energy(p);
  }

  ANY_DEVICE_INLINE double internal_energy(PressureEntropy theta) const {
    return ideal_internal_energy(theta.p());
  }

  ANY_DEVICE_INLINE double internal_energy(RhoP rhoP) const {
    return ideal_internal_energy(rhoP.p());
  }

  ANY_DEVICE_INLINE double internal_energy(RhoEntropy rhoK) const {
    double p = pressure(rhoK);
    return ideal_internal_energy(p);
  }

  /// -- Speed of sound.
  ANY_DEVICE_INLINE double sound_speed(const euler_var_t &u) const {
    return sound_speed(RhoP{u[0], pressure(u)});
  }

  ANY_DEVICE_INLINE double sound_speed(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    return zisa::sqrt(gamma() * p / rho);
  }

  ANY_DEVICE_INLINE double sound_speed(const RhoE &rhoE) const {
    return sound_speed(rhoP(rhoE));
  }

  // -- Temperature
  ANY_DEVICE_INLINE double temperature(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    return p / (rho * specific_gas_constant());
  }

  ANY_DEVICE_INLINE double temperature(const euler_var_t &u) const {
    return temperature(RhoP{u[0], pressure(u)});
  }

  ANY_DEVICE_INLINE double temperature(RhoE rhoE) const {
    return temperature(RhoP{rhoE.rho(), pressure(rhoE)});
  }

  euler_full_xvars_t full_extra_variables(const RhoE &rhoE) const {
    euler_full_xvars_t ret;

    auto [rho, E] = rhoE;
    auto p = pressure(rhoE);
    auto a = sound_speed(RhoP{rho, p});
    auto h = enthalpy(RhoP{rho, p});
    auto s = entropy(RhoP{rho, p});
    auto T = temperature(RhoP{rho, p});

    ret.rho = rho;
    ret.E = E;
    ret.a = a;
    ret.h = h;
    ret.s = s;
    ret.p = p;
    ret.T = T;

    return ret;
  }

  std::string str() const {
    return string_format(
        "Ideal gas EOS [%.3e, %.3e]", gamma(), specific_gas_constant());
  }

private:
  double gamma_;
  double specific_gas_constant_;
};

void save(HDF5Writer &writer, const IdealGasEOS &eos);

} // namespace zisa

#endif /* end of include guard */
