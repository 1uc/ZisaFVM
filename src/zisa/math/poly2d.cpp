/*
 *
 */

#include <zisa/math/poly2d.hpp>

namespace zisa {

int poly_degree(int_t n_coeffs) {

  // c = (d + 1)*(d +2) / 2
  // 0 = d**2 + 3*d + 2 - 2c
  // c_{1,2} = -3 +- sqrt(-4*(2-2c)) / 2

  LOG_ERR_IF(n_coeffs == 0, "Too few coefficients.");

  if (n_coeffs == 1) {
    return 0;
  }

  auto d1 = (-3.0 + zisa::sqrt(4.0 * (2.0 * n_coeffs - 2.0))) / 2.0;
  auto d2 = (-3.0 - zisa::sqrt(4.0 * (2.0 * n_coeffs - 2.0))) / 2.0;

  auto d = int(zisa::max(d1, d2));

  if (poly_dof(d + 1) <= n_coeffs) {
    return d + 1;
  } else if (poly_dof(d) <= n_coeffs) {
    return d;
  } else if (poly_dof(d - 1) <= n_coeffs) {
    return d - 1;
  } else {
    LOG_ERR("Broken logic.");
  }
}
} // namespace zisa
