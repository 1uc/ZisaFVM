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
  double alpha = PolytropeGravity::alpha();

  // Potential division by zero coming up:
  //     lim_{r->0} sin(r)/r = 1
  // The problem only exists for `r == 0`. We therefore add a small `eps`
  // which will almost always be truncated, except for the case of `r == 0`.
  double chi_eff = alpha * (chi + eps);
  return -2.0 * K * rhoC * zisa::sin(chi_eff) / chi_eff;
}

ANY_DEVICE_INLINE double PolytropeGravity::dphi_dx(double chi) const {
  double alpha = PolytropeGravity::alpha();
  double chi_eff = alpha * (chi + eps);

  double dphi = (zisa::cos(chi_eff) - zisa::sin(chi_eff) / chi_eff) / chi_eff;
  return -2.0 * K * rhoC * dphi * alpha;
}

ANY_DEVICE_INLINE double PolytropeGravity::alpha() const {
  return zisa::sqrt(2 * zisa::pi * G / K);
}

ANY_DEVICE_INLINE double PolytropeGravityRadial::alpha() const {
  return this->gravity.alpha();
}

// ---  SphericalGravity  -----------------------------------------------
ANY_DEVICE_INLINE double SphericalGravity::phi(double r) const {
  int_t i = index(r);

  double alpha = (r - radii(i)) / dr;
  return (1 - alpha) * phi_points[i] + alpha * phi_points[i + 1];
}

ANY_DEVICE_INLINE double SphericalGravity::dphi_dx(double r) const {
  int_t i = index(r);
  return (phi_points[i + 1] - phi_points[i]) / dr;
}

ANY_DEVICE_INLINE int_t SphericalGravity::index(double r) const {
  return int_t((r - domain[0]) / dr);
}

ANY_DEVICE_INLINE double SphericalGravity::radii(int_t i) const {
  return domain[0] + double(i) * dr;
}

} // namespace zisa
#endif /* end of include guard */
