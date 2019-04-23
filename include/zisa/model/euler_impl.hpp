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
template <class EOS, class Gravity>
Euler<EOS, Gravity>::Euler(const EOS &eos, const Gravity &gravity)
    : eos(eos), gravity(gravity) {}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE double
Euler<EOS, Gravity>::max_eigen_value(const euler_var_t &u) const {
  double p = eos.pressure(u);
  return zisa::sqrt((u(1) * u(1) + u(2) * u(2)) / (u(0) * u(0)))
         + eos.sound_speed(RhoP{u[0], p});
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE euler_var_t
Euler<EOS, Gravity>::flux(const euler_var_t &u) const {
  double p = eos.pressure(u);
  return flux(u, p);
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE euler_var_t Euler<EOS, Gravity>::flux(const euler_var_t &u,
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

template <class EOS, class Gravity>
ANY_DEVICE_INLINE double
Euler<EOS, Gravity>::energy(double rho, double v1, double v2, double p) const {
  return eos.internal_energy(RhoP{rho, p}) + 0.5 * rho * (v1 * v1 + v2 * v2);
}

template <class EOS, class Gravity>
void save(HDF5Writer &writer, const Euler<EOS, Gravity> &euler) {
  writer.open_group("model");
  writer.write_string("Euler", "name");
  save(writer, euler.eos);
  save(writer, euler.gravity);
  writer.close_group();
}

template <class EOS, class Gravity>
Euler<EOS, Gravity> Euler<EOS, Gravity>::load(HDF5Reader &reader) {
  reader.open_group("model");
  auto euler = Euler<EOS, Gravity>{EOS::load(reader), Gravity::load(reader)};
  reader.close_group();

  return euler;
}

template <class EOS, class Gravity>
std::string Euler<EOS, Gravity>::str() const {
  std::stringstream ss;

  ss << "Euler equations:\n"
     << indent_block(1, eos.str()) << "\n"
     << indent_block(1, gravity.str());

  return ss.str();
}

ANY_DEVICE_INLINE bool notplausible(const euler_var_t &u) {
  if (u(0) <= 0 || u(4) <= 0 || !isreal(u)) {
    return true;
  }

  return false;
}

ANY_DEVICE_INLINE bool isplausible(const euler_var_t &u) {
  if (u(0) <= 0 || u(4) <= 0 || !isreal(u)) {
    return false;
  }

  return true;
}

} // namespace zisa

#endif /* end of include guard: EULER_IMPL_CUH_5WLETFDJ */
