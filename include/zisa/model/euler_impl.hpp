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

namespace zisa {
template <class EOS, class Gravity>
Euler<EOS, Gravity>::Euler(double omega, const EOS &eos, const Gravity &gravity)
    : eos(eos), gravity(gravity), omega_(omega) {}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE double Euler<EOS, Gravity>::omega(void) const {
  return omega_;
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE double
Euler<EOS, Gravity>::max_eigen_value(const euler_var_t &u) const {
  double p = eos.pressure(u);
  return zisa::sqrt((u(1) * u(1) + u(2) * u(2)) / (u(0) * u(0)))
         + eos.sound_speed__rho_p(u[0], p);
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
ANY_DEVICE_INLINE euler_var_t Euler<EOS, Gravity>::physical_source(
    const euler_var_t &u, const XY &x, double /* t */) const {

  euler_var_t ps;
  double om = omega();

  double dphi_dx = gravity.dphi_dx(x, 0);
  double dphi_dy = gravity.dphi_dx(x, 1);

  // The coriolis force is
  //     coriolis = 2*(0, 0, omega) x (u(1), u(2), u(3))

  // density
  ps(0) = 0.0;

  // momentum
  ps(1) = -u(0) * dphi_dx + 2 * om * u(2);
  ps(2) = -u(0) * dphi_dy - 2 * om * u(1);
  ps(3) = 0.0;

  // energy
  ps(4) = -(u(1) * dphi_dx + u(2) * dphi_dy);

  return ps;
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE XY Euler<EOS, Gravity>::velocity(const euler_var_t &u) const {
  XY v;
  v[0] = u[1] / u[0];
  v[1] = u[2] / u[0];
  return v;
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE double
Euler<EOS, Gravity>::energy(double rho, double v1, double v2, double p) const {
  return eos.internal_energy__rho_p(rho, p) + 0.5 * rho * (v1 * v1 + v2 * v2);
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE euler_var_t Euler<EOS, Gravity>::natural_variables(
    double rho, double v1, double v2, double p) const {
  euler_var_t ret;

  ret(0) = rho;
  ret(1) = rho * v1;
  ret(2) = rho * v2;
  ret(3) = 0.0;
  ret(4) = energy(rho, v1, v2, p);

  return ret;
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE euler_var_t Euler<EOS, Gravity>::natural_variables(
    double rho, const XY &v, double p) const {
  return natural_variables(rho, v[0], v[1], p);
}

template <class EOS, class Gravity>
ANY_DEVICE_INLINE euler_var_t
Euler<EOS, Gravity>::primitive_variables(const euler_var_t &u_in) const {
  euler_var_t u_out;

  double p = eos.pressure(u_in);

  // rho (unchanged)
  u_out(0) = u_in(0);

  // momentum --> velocity
  u_out(1) = u_in(1) / u_in(0);
  u_out(2) = u_in(2) / u_in(0);
  u_out(3) = 0.0;

  // energy --> pressure
  u_out(4) = p;

  return u_out;
}

namespace detail {

inline void save_parameters(const IdealGasEOS &eos, HDF5Writer &writer) {
  writer.write_string("IdealGasEOS", "eos");
  writer.write_scalar(eos.gamma(), "gamma");
  writer.write_scalar(eos.specific_gas_constant(), "specific_gas_constant");
}

} // namespace detail

template <class EOS, class Gravity>
void Euler<EOS, Gravity>::save_parameters(HDF5Writer &writer) const {
  writer.write_string("Euler", "model");
  writer.write_scalar(omega(), "omega");

  detail::save_parameters(eos, writer);
}

// template <class EOS, class Gravity>
// void Euler<EOS, Gravity>::load_state(AllVariables &u, HDF5Reader &reader) {
//   int n_advected_variables = reader.read_scalar<int>("n_advected_variables");
//   assert(u.advected_variables.shape[2] == n_advected_variables);

//   u.conserved_variables.split_load(reader, all_labels<cvars_t>());
//   u.advected_variables.split_load(
//       reader, numbered_labels("mq%d", n_advected_variables));
// }

template <class EOS, class Gravity>
std::string Euler<EOS, Gravity>::str(int indent) const {
  std::stringstream ss;
  std::string space(2 * indent, ' ');
  std::string space2(2 * (indent + 1), ' ');

  ss << space << "Euler equations:\n";
  ss << space2 << string_format("omega   = %.3e \n", omega());

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
