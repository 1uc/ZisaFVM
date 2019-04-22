#ifndef ZISA_JANKA_EOS_2398ESQ_HPP
#define ZISA_JANKA_EOS_2398ESQ_HPP

#include <array>

#include <zisa/cli/input_parameters.hpp>
#include <zisa/math/brent.hpp>
#include <zisa/math/newton.hpp>
#include <zisa/model/equation_of_state.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

struct JankaEOSParams {
  double rho_bounce;
  std::array<double, 2> gamma;
  double gamma_thermal;
  double E1;

  std::string str() const;
};
JankaEOSParams make_default_janka_eos_params();

void save(HDF5Writer &writer, const JankaEOSParams &params);

class JankaEOS : public EquationOfState {
private:
  using cvars_t = euler_var_t;
  using xvars_t = cvars_t::xvars_t;

public:
  JankaEOS() = default;
  JankaEOS(double rho_bounce,
           const std::array<double, 2> &gamma,
           double gamma_thermal,
           double E1)
      : rho_bounce(rho_bounce), gamma(gamma), gamma_thermal(gamma_thermal) {

    E[0] = E1;
    E[1] = (gamma[0] - 1.0) / (gamma[1] - 1.0) * E1
           * pow(rho_bounce, gamma[0] - gamma[1]);

    E[2] = (gamma[1] - gamma[0]) / (gamma[1] - 1.0) * E1
           * pow(rho_bounce, gamma[0] - 1.0);

    K[0] = (gamma[0] - 1.0) * E[0];
    K[1] = (gamma[1] - 1.0) * E[1];

    h_rho_bounce = enthalpy_fixed<0>(rho_bounce);
  }

  JankaEOS(const JankaEOSParams &params)
      : JankaEOS(
            params.rho_bounce, params.gamma, params.gamma_thermal, params.E1) {}

  JankaEOSParams params() const {
    return {rho_bounce, gamma, gamma_thermal, E[0]};
  }

  using EquationOfState::internal_energy;
  using EquationOfState::kinetic_energy;

  // --- Ideal Gas EOS calls ---------------------------------------------------
  double ideal_gas_enthalpy(double K, double rho, double gamma) const {
    return gamma / (gamma - 1.0) * K * pow(rho, gamma - 1.0);
  }

  /// Returns the ideal gas pressure and it's derivative w.r.t. rho.
  double ideal_gas_pressure(double K, double rho, double gamma) const {
    return K * zisa::pow(rho, gamma);
  }

  // --- gamma  ----------------------------------------------------------------
  double polytropic_gamma(double rho) const {
    return rho <= rho_bounce ? gamma[0] : gamma[1];
  }

  // --- Density ---------------------------------------------------------------
  double ideal_gas_density(double gamma, const EnthalpyEntropy &theta) const {
    const auto &[h, K] = theta;
    return zisa::pow((gamma - 1) / gamma * h / K, 1.0 / (gamma - 1));
  }

  template <int REGIME>
  double rho_fixed(EnthalpyEntropy theta) const {
    const auto &[h, K] = theta;

    double K_p = this->K[REGIME];
    double K_th = K - K_p;

    auto f = [this, h = h, K_th](double rho) {
      double h_p = enthalpy_fixed<REGIME>(rho);
      double h_th = ideal_gas_enthalpy(K_th, rho, gamma_thermal);

      return h - (h_p + h_th);
    };

    double rho_a = ideal_gas_density(gamma[REGIME], theta);
    double rho_b = 2.0 * rho_a; // fixme this.

    double inf = std::numeric_limits<double>::infinity();
    auto hard_bounds = (REGIME == 0 ? std::array<double, 2>{0.0, rho_bounce}
                                    : std::array<double, 2>{rho_bounce, inf});

    double atol = 1e-10 * zisa::avg(rho_a, rho_b);
    int max_iter = 100;

    std::tie(rho_a, rho_b)
        = find_bracket(f, {rho_a, rho_b}, hard_bounds, max_iter);

    return brent(f, rho_a, rho_b, atol, max_iter);
  }

  ANY_DEVICE_INLINE double rho(EnthalpyEntropy theta) const {
    const auto &[h, K] = theta;

    if (h <= h_rho_bounce) {
      // definitely rho < rho_b
      return this->rho_fixed<0>(theta);
    } else {
      // definitely rho > rho_b
      return this->rho_fixed<1>(theta);
    }
  }

  // --- Pressure  -------------------------------------------------------------
  ANY_DEVICE_INLINE double polytropic_pressure(double rho) const {
    int_t i = rho <= rho_bounce ? 0 : 1;
    return K[i] * zisa::pow(rho, gamma[i]);
  }

  ANY_DEVICE_INLINE double thermal_pressure(double thermal_energy) const {
    return (gamma_thermal - 1.0) * thermal_energy;
  }

  ANY_DEVICE_INLINE double thermal_pressure(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    return p - polytropic_pressure(rho);
  }

  ANY_DEVICE_INLINE double total_pressure(double rho, double E_th) const {
    return polytropic_pressure(rho) + thermal_pressure(E_th);
  }

  ANY_DEVICE_INLINE double pressure(RhoE rhoE) const {
    const auto &[rho, E] = rhoE;
    double E_th = thermal_energy(rhoE);
    return total_pressure(rho, E_th);
  }

  ANY_DEVICE_INLINE double pressure(RhoEntropy rhoK) const {
    const auto &[rho, K] = rhoK;

    double K_thermal = K - polytropic_entropy(rho);
    double p_thermal = ideal_gas_pressure(K_thermal, rho, gamma_thermal);

    return polytropic_pressure(rho) + p_thermal;
  }

  ANY_DEVICE_INLINE double pressure(EnthalpyEntropy theta) const {
    double rho = this->rho(theta);
    return pressure(RhoEntropy{rho, theta.s()});
  }

  ANY_DEVICE_INLINE double pressure(const cvars_t &u) const {
    return pressure(RhoE{u[0], internal_energy(u)});
  }

  ANY_DEVICE_INLINE std::pair<double, double>
  split_pressure(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    double p_polytrope = polytropic_pressure(rho);
    return {p_polytrope, p - p_polytrope};
  }

  // --- Energy ----------------------------------------------------------------
  ANY_DEVICE_INLINE double thermal_energy(double p_thermal) const {
    return p_thermal / (gamma_thermal - 1.0);
  }

  ANY_DEVICE_INLINE double thermal_energy(RhoE rhoE) const {
    const auto &[rho, E] = rhoE;
    return E - polytropic_energy(rho);
  }

  ANY_DEVICE_INLINE double polytropic_energy(double rho) const {
    if (rho <= rho_bounce) {
      return E[0] * pow(rho, gamma[0]);
    } else {
      return E[1] * pow(rho, gamma[1]) + E[2] * rho;
    }
  }

  ANY_DEVICE_INLINE double internal_energy(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;

    double p_thermal = thermal_pressure(rhoP);
    return polytropic_energy(rho) + thermal_energy(p_thermal);
  }

  ANY_DEVICE_INLINE double internal_energy(RhoE rhoE) const { return rhoE.E(); }
  ANY_DEVICE_INLINE double internal_energy(RhoEntropy rhoK) const {
    const auto &[rho, K] = rhoK;

    return internal_energy(RhoP{rho, pressure(rhoK)});
  }

  // --- Enthalpy --------------------------------------------------------------
  template <int REGIME>
  double enthalpy_fixed(double rho) const {
    double h_p = ideal_gas_enthalpy(K[REGIME], rho, gamma[REGIME]);
    return h_p + (REGIME == 0 ? 0.0 : E[2]);
  }

  double polytropic_enthalpy(double rho) const {
    if (rho <= rho_bounce) {
      return enthalpy_fixed<0>(rho);
    } else {
      return enthalpy_fixed<1>(rho);
    }
  }

  ANY_DEVICE_INLINE double enthalpy(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    double E = internal_energy(rhoP);
    return (E + p) / rho;
  }

  ANY_DEVICE_INLINE double enthalpy(const RhoEntropy &rhoK) const {
    return enthalpy(rhoP(rhoK));
  }

  // --- Entropy ---------------------------------------------------------------
  ANY_DEVICE_INLINE double polytropic_entropy(double rho) const {
    return rho <= rho_bounce ? K[0] : K[1];
  }

  ANY_DEVICE_INLINE double entropy(const RhoE &rhoE) const {
    return entropy(rhoP(rhoE));
  }

  ANY_DEVICE_INLINE double entropy(const RhoP &rhoP) const {
    const auto &[rho, p] = rhoP;
    const auto &[p_polytropic, p_thermal] = split_pressure(rhoP);

    double K_polytropic = polytropic_entropy(rho);
    double K_thermal = p_thermal * zisa::pow(rho, -gamma_thermal);

    return K_polytropic + K_thermal;
  }

  // --- speed of sound --------------------------------------------------------
  ANY_DEVICE_INLINE double sound_speed(const RhoE &rhoE) const {
    return sound_speed(rhoP(rhoE));
  }

  ANY_DEVICE_INLINE double sound_speed(const RhoP &rhoP) const {
    double rho = rhoP.rho();
    double K = entropy(rhoP);
    double K_p = polytropic_entropy(rho);
    double K_th = K - K_p;

    double gamma_p = polytropic_gamma(rho);

    double cs2 = gamma_p * K_p * zisa::pow(rho, gamma_p - 1.0)
                 + gamma_thermal * K_th * zisa::pow(rho, gamma_thermal - 1.0);

    return zisa::sqrt(cs2);
  }

  // -- 2 -> 2 Conversions -----------------------------------------------------
  ANY_DEVICE_INLINE RhoP rhoP(const cvars_t &u) const {
    return RhoP{u[0], pressure(u)};
  }

  ANY_DEVICE_INLINE RhoP rhoP(const RhoE &rhoE) const {
    return RhoP{rhoE.rho(), pressure(rhoE)};
  }

  ANY_DEVICE_INLINE RhoP rhoP(const RhoEntropy &rhoK) const {
    return RhoP{rhoK.rho(), pressure(rhoK)};
  }

  ANY_DEVICE_INLINE RhoP rhoP(const EnthalpyEntropy &theta) const {
    double rho = this->rho(theta);
    return RhoP{rho, pressure(RhoEntropy{rho, theta.s()})};
  }

  ANY_DEVICE_INLINE RhoE rhoE(const RhoP &rhoP) const {
    return RhoE{rhoP.rho(), internal_energy(rhoP)};
  }

  ANY_DEVICE_INLINE RhoE rhoE(const RhoEntropy &rhoK) const {
    return RhoE{rhoK.rho(), internal_energy(rhoK)};
  }

  ANY_DEVICE_INLINE RhoE rhoE(const EnthalpyEntropy &theta) const {
    double rho = this->rho(theta);
    return {rho, internal_energy(RhoEntropy{rho, theta.s()})};
  }

  ANY_DEVICE_INLINE RhoE rhoE(const cvars_t &u) const {
    return RhoE{u[0], internal_energy(u)};
  }

  ANY_DEVICE_INLINE EnthalpyEntropy enthalpy_entropy(const RhoP &rhoP) const {
    return {enthalpy(rhoP), entropy(rhoP)};
  }

  ANY_DEVICE_INLINE EnthalpyEntropy enthalpy_entropy(const RhoE &rhoE) const {
    return enthalpy_entropy(rhoP(rhoE));
  }

  ANY_DEVICE_INLINE EnthalpyEntropy
  enthalpy_entropy(const RhoEntropy &rhoK) const {
    return {enthalpy(rhoP(rhoK)), rhoK.s()};
  }

  ANY_DEVICE_INLINE RhoEntropy rhoK(const RhoP &rhoP) const {
    return {rhoP.rho(), entropy(rhoP)};
  }

  ANY_DEVICE_INLINE RhoEntropy rhoK(const RhoE &rhoE) const {
    return {rhoE.rho(), entropy(rhoE)};
  }

  ANY_DEVICE_INLINE RhoEntropy rhoK(const EnthalpyEntropy &theta) const {
    return {rho(theta), theta.s()};
  }

  // --- cvars, xvars, etc. ----------------------------------------------------
  ANY_DEVICE_INLINE cvars_t cvars(RhoE rhoE) const {
    return {rhoE.rho(), 0.0, 0.0, 0.0, rhoE.E()};
  }

  ANY_DEVICE_INLINE xvars_t xvars(const cvars_t &u) const {
    return {pressure(u), sound_speed(rhoE(u))};
  }

  std::string str() const;

private:
  double rho_bounce;
  double h_rho_bounce;
  std::array<double, 2> gamma;
  double gamma_thermal;
  std::array<double, 3> E;
  std::array<double, 2> K;
};

void save(HDF5Writer &writer, const JankaEOS &eos);

JankaEOS make_janka_eos(const JankaEOSParams &params);

JankaEOS make_default_janka_eos();

}
#endif // ZISA_JANKA_EOS_HPP
