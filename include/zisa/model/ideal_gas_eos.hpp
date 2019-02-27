/* Ideal Gas equation of state.
 */

#ifndef IDEAL_GAS_EOS_H_SNFCI
#define IDEAL_GAS_EOS_H_SNFCI

#include <zisa/model/equation_of_state_base.hpp>

namespace zisa {

class IdealGasEOS : public EquationOfState {
private:
  using cvars_t = euler_var_t;
  using xvars_t = cvars_t::xvars_t;

public:
  IdealGasEOS() = default;
  IdealGasEOS(double gamma, double specific_gas_constant)
      : gamma_(gamma), specific_gas_constant_(specific_gas_constant) {}

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
    auto [rho, E] = rhoE;
    return cvars_t{rho, 0.0, 0.0, 0.0, E};
  }

  ANY_DEVICE_INLINE EnthalpyEntropy enthalpy_entropy(RhoE rhoE) const {

    double p = pressure__rho_E(rhoE.rho(), rhoE.E());
    double h = enthalpy__rho_p(rhoE.rho(), p);
    double K = K__rho_p(rhoE.rho(), p);

    return EnthalpyEntropy{h, K};
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

  ANY_DEVICE_INLINE RhoT rhoT(RhoE rhoE) const {

    auto rho = rhoE.rho();
    auto T = this->temperature(rhoE);

    return RhoT{rho, T};
  }

  ANY_DEVICE_INLINE double entropy(RhoE rhoE) const { return K(rhoE); }

  ANY_DEVICE_INLINE double entropy__rho_p(double rho, double p) const {
    return cV() * zisa::log(p / zisa::pow(rho, gamma()));
  }

  ANY_DEVICE_INLINE double enthalpy__rho_p(double rho, double p) const {
    return gamma() / (gamma() - 1.0) * p / rho;
  }

  ANY_DEVICE_INLINE double enthalpy(const RhoE &rhoE) const {
    return enthalpy__rho_p(rhoE.rho(), pressure(rhoE));
  }

  ANY_DEVICE_INLINE double K__rho_p(double rho, double p) const {
    return p / zisa::pow(rho, gamma());
  }

  ANY_DEVICE_INLINE double K(const RhoE &rhoE) const {
    return K__rho_p(rhoE.rho(), pressure(rhoE));
  }

  ANY_DEVICE_INLINE double rho__p_T(double p, double T) const {
    return p / (T * specific_gas_constant());
  }

  ANY_DEVICE_INLINE double rho__h_K(double h, double K) const {
    double gamma = IdealGasEOS::gamma();
    double base = 1.0 / K * (gamma - 1.0) / gamma * h;
    double exponent = 1.0 / (gamma - 1.0);
    return zisa::pow(base, exponent);
  }

  ANY_DEVICE_INLINE double rho__p_K(double p, double K) const {
    return zisa::pow(p / K, 1.0 / gamma());
  }

  ANY_DEVICE_INLINE double rho__p_s(double p, double s) const {
    double gamma = IdealGasEOS::gamma();
    double cV = IdealGasEOS::cV();

    return zisa::pow(p, 1.0 / gamma) * zisa::exp(-s / (cV * gamma));
  }

  ANY_DEVICE_INLINE double rho(EnthalpyEntropy theta) const {
    return rho__h_K(theta.h(), theta.K());
  }

  ANY_DEVICE_INLINE double rho(PressureEntropy theta) const {
    return rho__p_K(theta.p(), theta.s());
  }

  /// Pressure given density and internal energy, using ideal gas law.
  ANY_DEVICE_INLINE double pressure__rho_E(double, double E) const {
    return E * (gamma() - 1.0);
  }

  ANY_DEVICE_INLINE double pressure__rho_K(double rho, double K) const {
    return K * zisa::pow(rho, gamma());
  }

  ANY_DEVICE_INLINE double pressure__rho_T(double rho, double T) const {
    double r_gas = specific_gas_constant();
    return r_gas * rho * T;
  }

  ANY_DEVICE_INLINE double pressure__rho_s(double rho, double s) const {
    // s = cV() * zisa::log(p / zisa::pow(rho, gamma()));
    return zisa::exp(s / cV()) * zisa::pow(rho, gamma());
  }

  ANY_DEVICE_INLINE double pressure(RhoE rhoE) const {
    return pressure__rho_E(rhoE.rho(), rhoE.E());
  }

  ANY_DEVICE_INLINE double pressure(const euler_var_t &u) const {
    double E_kin = kinetic_energy(u);
    return pressure(RhoE{u[0], u[4] - E_kin});
  }

  /// Kinetic energy.
  ANY_DEVICE_INLINE double kinetic_energy(const euler_var_t &u) const {
    return zisa::kinetic_energy(u);
  }

  /// Internal energy.
  ANY_DEVICE_INLINE double internal_energy__p(double p) const {
    return p / (gamma() - 1.0);
  }

  ANY_DEVICE_INLINE double internal_energy__rho_p(double, double p) const {
    return internal_energy__p(p);
  }

  ANY_DEVICE_INLINE double internal_energy__h_K(double h, double K) const {
    double rho = rho__h_K(h, K);
    double p = pressure__rho_K(rho, K);
    return internal_energy__rho_p(rho, p);
  }

  ANY_DEVICE_INLINE double internal_energy__rho_K(double rho, double K) const {
    double p = pressure__rho_K(rho, K);
    return internal_energy__rho_p(rho, p);
  }

  ANY_DEVICE_INLINE double internal_energy(const euler_var_t &u) const {
    return zisa::internal_energy(u);
  }

  ANY_DEVICE_INLINE double internal_energy(EnthalpyEntropy theta) const {
    return internal_energy__h_K(theta.h(), theta.K());
  }

  ANY_DEVICE_INLINE double internal_energy(PressureEntropy theta) const {
    return internal_energy__p(theta.p());
  }

  ANY_DEVICE_INLINE double internal_energy(RhoP rhoP) const {
    return internal_energy__p(rhoP.p());
  }

  /// Specific internal energy.
  ANY_DEVICE_INLINE double specific_internal_energy__rho_E(double rho,
                                                           double E) const {
    return E / rho;
  }

  /// -- Speed of sound.
  ANY_DEVICE_INLINE double sound_speed(const euler_var_t &u) const {
    return sound_speed__rho_p(u[0], pressure(u));
  }

  ANY_DEVICE_INLINE double sound_speed__rho_p(double rho, double p) const {
    return zisa::sqrt(gamma() * p / rho);
  }

  // -- Temperature
  ANY_DEVICE_INLINE double temperature__rho_p(double rho, double p) const {
    return p / (rho * specific_gas_constant());
  }

  ANY_DEVICE_INLINE double temperature(const euler_var_t &u) const {
    return temperature__rho_p(u[0], pressure(u));
  }

  ANY_DEVICE_INLINE double temperature(RhoE rhoE) const {
    return temperature__rho_p(rhoE.rho(), pressure(rhoE));
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
