/* Implementation details, fake a .cu file (there is a .cu file, but these are
 * all inline functions).
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2014-12-21
 */
#ifndef EULER_IMPL_CUH_5WLETFDJ
#define EULER_IMPL_CUH_5WLETFDJ

#include "euler_decl.hpp"

#include <zisa/math/basic_functions.hpp>
#include <zisa/math/isreal.hpp>
#include <zisa/model/equation_of_state.hpp>
#include <zisa/utils/indent_block.hpp>

namespace zisa {
ANY_DEVICE_INLINE double Euler::max_eigen_value(const euler_var_t &u,
                                                double cs) const {

  double v2 = (u(1) * u(1) + u(2) * u(2) + u(3) * u(3)) / (u(0) * u(0));

  return std::sqrt(v2) + cs;
}

ANY_DEVICE_INLINE euler_var_t Euler::flux(const euler_var_t &u,
                                          double p) const {
  euler_var_t pf;

  double v = u(1) / u(0);

  pf(0) = u(1);
  pf(1) = v * u(1) + p;
  pf(2) = v * u(2);
  pf(3) = v * u(3);
  pf(4) = v * (u(4) + p);

  return pf;
}

template <class EOS>
ANY_DEVICE_INLINE double
total_energy(const EOS &eos, double rho, double v1, double v2, double p) {
  return eos.internal_energy(RhoP{rho, p}) + 0.5 * rho * (v1 * v1 + v2 * v2);
}

ANY_DEVICE_INLINE bool notplausible(const euler_var_t &u) {
  return u(0) <= 0 || u(4) <= 0 || !isreal(u);
}

ANY_DEVICE_INLINE bool isplausible(const euler_var_t &u) {
  return !notplausible(u);
}

} // namespace zisa

#endif /* end of include guard: EULER_IMPL_CUH_5WLETFDJ */
