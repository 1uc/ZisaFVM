/* Callbacks for the Euler equations.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2014-11-12
 */
#ifndef EULER_CUH_OTPMLDSW
#define EULER_CUH_OTPMLDSW

#include <zisa/config.hpp>

#include <zisa/io/hdf5_writer.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/model/equation_of_state.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/models.hpp>

namespace zisa {
/// Implements a `Model` for the Euler equations.
template <class EOS, class Gravity>
class Euler {
public:
  using eos_t = EOS;
  using gravity_t = Gravity;
  using cvars_t = euler_var_t;

public:
  const eos_t eos;
  const gravity_t gravity;

public:
  Euler() = default;
  Euler(const EOS &eos, const Gravity &gravity);

  /// Largest eigenvalue of the derivative of the flux.
  ANY_DEVICE_INLINE double max_eigen_value(const euler_var_t &u) const;

  /// Physical flux normal to the surface.
  /** The convention is that the first component is normal to the surface
   *  through which the flux is computed.
   */
  ANY_DEVICE_INLINE euler_var_t flux(const euler_var_t &u) const;
  ANY_DEVICE_INLINE euler_var_t flux(const euler_var_t &u, double p) const;

  /// Total energy.
  ANY_DEVICE_INLINE double
  energy(double rho, double v1, double v2, double p) const;

  /// Kinetic energy of the fluid parcel.
  ANY_DEVICE_INLINE double kinetic_energy(const euler_var_t &u) const;

  /// Self-documenting string.
  std::string str() const;
};

/// Write the parameters of this model to disk.
template<class EOS, class Gravity>
void save(HDF5Writer &writer, const Euler<EOS, Gravity> &euler);


/// Are the value obviously unphysical?
/** The variable `u` is considered implausible if any component is not even
 *  a real number or if the density, (pressure, GPU only) or energy is
 *  non-positive.
 *
 *  @param u variable to check
 *  @return true if clearly unphysical.
 */
ANY_DEVICE_INLINE bool notplausible(const euler_var_t &u);

/// Could these values be physical (at first glace)?
/** The variable `u` is considered plausible if all its components are real
 *  numbers and furthermore, density, (pressure, GPU only) and energy are
 *  (strictly) positive.
 *
 *  @param u variable to check
 *  @return true if clearly unphysical.
 */
ANY_DEVICE_INLINE bool isplausible(const euler_var_t &u);

template <class Model>
class euler_like {
public:
  static constexpr bool value = false;
};

template <class EOS, class Gravity>
class euler_like<Euler<EOS, Gravity>> {
public:
  static constexpr bool value = true;
};

} // namespace zisa

#endif /* end of include guard: EULER_H_OTPMLDSW */
