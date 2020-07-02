/* The famous HLLC numerical flux.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2015-02-18
 */
#ifndef HLLC_H_WOCXPIIT
#define HLLC_H_WOCXPIIT

#include <zisa/model/equation_of_state.hpp>
#include <zisa/model/euler.hpp>

#if ZISA_HAS_HELMHOLTZ_EOS == 1
#include <zisa/model/helmholtz_eos.hpp>
#endif

namespace zisa {

/// The Roe average of two states.
class RoeAverage {
public:
  /// Construct a RoeAverage w.r.t. the densities `rho_L` and `rho_R`.
  ANY_DEVICE_INLINE RoeAverage(double rho_L, double rho_R)
      : roe_ratio(zisa::sqrt(rho_R / rho_L)) {}

  /// Compute the Roe average.
  ANY_DEVICE_INLINE double operator()(double qL, double qR) {
    return (qL + qR * roe_ratio) / (1.0 + roe_ratio);
  }

private:
  double roe_ratio;
};

template <class Model>
class hllc_speeds;

template <class Gravity>
class hllc_speeds<Euler<IdealGasEOS, Gravity>> {
private:
  using euler_t = Euler<IdealGasEOS, Gravity>;
  using cvars_t = euler_var_t;
  using xvars_t = euler_var_t::xvars_t;

public:
  static std::tuple<double, double, double> speeds(const euler_t &euler,
                                                   const cvars_t &uL,
                                                   const xvars_t &xvarL,
                                                   const cvars_t &uR,
                                                   const xvars_t &xvarR) {

    RoeAverage roe_average(uL(0), uR(0));

    auto pL = xvarL.p;
    auto aL = xvarL.a;

    auto pR = xvarR.p;
    auto aR = xvarR.a;

    double vL = uL(1) / uL(0);
    double vR = uR(1) / uR(0);
    double v_tilda = roe_average(vL, vR);

    double HL = (uL(4) + pL) / uL(0);
    double HR = (uR(4) + pR) / uR(0);
    double H_tilda = roe_average(HL, HR);

    double vroe_square
        = zisa::pow<2>(roe_average(uL(1) / uL(0), uR(1) / uR(0)))
          + zisa::pow<2>(roe_average(uL(2) / uL(0), uR(2) / uR(0)))
          + zisa::pow<2>(roe_average(uL(3) / uL(0), uR(3) / uR(0)));

    double a_tilda
        = zisa::sqrt((euler.eos.gamma() - 1.0) * (H_tilda - 0.5 * vroe_square));

    double sL = zisa::min(vL - aL, v_tilda - a_tilda);
    double sR = zisa::max(vR + aR, v_tilda + a_tilda);
    double s_star = (uR(1) * (sR - vR) - uL(1) * (sL - vL) + pL - pR)
                    / (uR(0) * (sR - vR) - uL(0) * (sL - vL));

    return {sL, s_star, sR};
  }
};

#if ZISA_HAS_HELMHOLTZ_EOS == 1

template <class Gravity>
class hllc_speeds<Euler<HelmholtzEOS, Gravity>> {
private:
  using euler_t = Euler<HelmholtzEOS, Gravity>;
  using cvars_t = euler_var_t;
  using xvars_t = euler_var_t::xvars_t;

public:
  static std::tuple<double, double, double> speeds(const euler_t &,
                                                   const cvars_t &uL,
                                                   const xvars_t &xvarL,
                                                   const cvars_t &uR,
                                                   const xvars_t &xvarR) {

    auto pL = xvarL.p;
    auto aL = xvarL.a;

    auto pR = xvarR.p;
    auto aR = xvarR.a;

    double vL = uL(1) / uL(0);
    double vR = uR(1) / uR(0);

    double sL = zisa::min(vL - aL, vR - aR);
    double sR = zisa::max(vL + aL, vR + aR);
    double s_star = (uR(1) * (sR - vR) - uL(1) * (sL - vL) + pL - pR)
                    / (uR(0) * (sR - vR) - uL(0) * (sL - vL));

    return {sL, s_star, sR};
  }
};

#endif

/// HLLC numerical flux with Einfeldt-Batten wavespeeds.
/** Reference: Batten, Wavespeed Estimates for the HLLC Riemann Solver, 1997
 */
template <class Model>
class HLLCBatten {
private:
public:
  using cvars_t = typename Model::cvars_t;
  using model_t = Model;
  using speeds_t = std::tuple<double, double, double>;

public:
  /// Compute the numerical flux.
  /*
   *  @param euler
   *  @param uL
   *    conserved variables 'inside' of the cell, w.r.t the normal
   *  @param uR
   *    conserved variables 'outside' of the cell, w.r.t. the normal.
   */
  ANY_DEVICE_INLINE static std::tuple<euler_var_t, speeds_t>
  flux(const model_t &euler, const cvars_t &uL, const cvars_t &uR) {

    const auto xvarL = euler.eos.xvars(uL);
    const auto xvarR = euler.eos.xvars(uR);

    const auto [sL, s_star, sR]
        = hllc_speeds<model_t>::speeds(euler, uL, xvarL, uR, xvarR);

    // HLLC flux
    const cvars_t &uK = (0.0 <= s_star ? uL : uR);
    const double pK = (0.0 <= s_star ? xvarL.p : xvarR.p);
    auto nf = euler.flux(uK, pK);

    if (sL < 0.0 && 0.0 <= sR) {
      double sK = (0.0 <= s_star ? sL : sR);
      double vK = (0.0 <= s_star ? uL(1) / uL(0) : uR(1) / uR(0));
      double cK = (sK - vK) / (sK - s_star);

      nf(0) += sK * (cK * uK(0) - uK(0));
      nf(1) += sK * (cK * uK(0) * s_star - uK(1));
      nf(2) += sK * (cK * uK(2) - uK(2));
      nf(3) += sK * (cK * uK(3) - uK(3));
      nf(4)
          += sK
             * (cK * (uK(4) + (s_star - vK) * (uK(0) * s_star + pK / (sK - vK)))
                - uK(4));
    }

    return {nf, std::tuple(sL, s_star, sR)};
  }

  ANY_DEVICE_INLINE static double tracer_flux(const cvars_t &uL,
                                              const cvars_t &uR,
                                              double mqL,
                                              double mqR,
                                              const speeds_t &speeds) {

    const auto [sL, s_star, sR] = speeds;
    double mqK = (0.0 <= s_star ? mqL : mqR);
    double vK = (0.0 <= s_star ? uL(1) / uL(0) : uR(1) / uR(0));
    double fK = mqK * vK;

    if (sL < 0.0 && 0.0 < sR) {
      double sK = (0.0 <= s_star ? sL : sR);
      double cK = (sK - vK) / (sK - s_star);

      return fK + sK * (cK * mqK - mqK);
    }

    return fK;
  }

  /// Self-documenting string.
  static std::string str() {
    return "HLLC with Batten wavespeeds (`HLLCBatten`)";
  }
};

} // namespace zisa

#endif /* end of include guard */
