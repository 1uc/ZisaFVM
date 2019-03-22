#ifndef GRAVITY_INC_H01YR
#define GRAVITY_INC_H01YR

#include "gravity.hpp"

#include <zisa/model/euler.hpp>

namespace zisa {
// --- ConstantGravity ---------------------------------------------
ANY_DEVICE_INLINE double ConstantGravity::phi(double chi) const {
  return gravity * chi;
}

ANY_DEVICE_INLINE double ConstantGravity::dphi_dx(double /* chi */) const {
  return gravity;
}

// --- PointMassGravity ---------------------------------------------
ANY_DEVICE_INLINE double PointMassGravity::phi(double chi) const {
  return GM / (X + chi);
}

ANY_DEVICE_INLINE double PointMassGravity::dphi_dx(double chi) const {
  return -GM / zisa::pow<2>(X + chi);
}

// ---  Polytrope  ---------------------------------------------------
inline PolytropeGravity::PolytropeGravity(double rhoC, double K, double G)
    : rhoC(rhoC), K(K), G(G) {}

ANY_DEVICE_INLINE double PolytropeGravity::phi(double chi) const {
  double alpha = PolytropeGravity::alpha(chi);

  // Potential division by zero coming up:
  //     lim_{r->0} sin(r)/r = 1
  // The problem only exists for `r == 0`. The correct limiting behaviour
  // is guaranteed. we therefore add a small `eps` which will almost always be
  // truncated, except for the case of `r == 0`.
  double chi_eff = alpha * (chi + eps);
  return -2.0 * K * rhoC * zisa::sin(chi_eff) / chi_eff;
}

ANY_DEVICE_INLINE double PolytropeGravity::dphi_dx(double chi) const {
  double alpha = PolytropeGravity::alpha(chi);
  double chi_eff = alpha * (chi + eps);

  double dphi = (zisa::cos(chi_eff) - zisa::sin(chi_eff) / chi_eff) / chi_eff;
  return -2.0 * K * rhoC * dphi * alpha;
}

ANY_DEVICE_INLINE double PolytropeGravity::alpha(double /* chi */) const {
  return zisa::sqrt(2 * zisa::pi * G / K);
}

ANY_DEVICE_INLINE double PolytropeGravityRadial::alpha(double chi) const {
  return this->gravity.alpha(chi);
}

// ---  PolytropeWithJump  ----------------------------------------------
inline PolytropeGravityWithJump::PolytropeGravityWithJump(
    double r_crit, double rhoC, double K_inner, double K_outer, double G)
    : r_crit(r_crit), inner(rhoC, K_inner, G), outer(1.0, K_outer, G) {

  double rhoC_outer = inner.phi(r_crit) / outer.phi(r_crit);
  outer = PolytropeGravity(rhoC_outer, K_outer, G);
}

ANY_DEVICE_INLINE double PolytropeGravityWithJump::phi(double chi) const {
  return (chi < r_crit) ? inner.phi(chi) : outer.phi(chi);
}

ANY_DEVICE_INLINE double PolytropeGravityWithJump::dphi_dx(double chi) const {
  return (chi < r_crit) ? inner.dphi_dx(chi) : outer.dphi_dx(chi);
}

ANY_DEVICE_INLINE double PolytropeGravityWithJump::alpha(double chi) const {
  return (chi < r_crit) ? inner.alpha(chi) : outer.alpha(chi);
}

ANY_DEVICE_INLINE double
PolytropeGravityWithJumpRadial::alpha(double chi) const {
  return this->gravity.alpha(chi);
}

} // namespace zisa
#endif /* end of include guard */
