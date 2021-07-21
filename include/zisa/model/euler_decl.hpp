// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

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
class Euler {
public:
  using cvars_t = euler_var_t;

public:
  Euler() = default;

  /// Largest eigenvalue of the derivative of the flux.
  ANY_DEVICE_INLINE double max_eigen_value(const euler_var_t &u,
                                           double cs) const;

  /// Physical flux normal to the surface.
  /** The convention is that the first component is normal to the surface
   *  through which the flux is computed.
   */
  ANY_DEVICE_INLINE euler_var_t flux(const euler_var_t &u, double p) const;

  /// Self-documenting string.
  std::string str() const;
};

template <class EOS>
ANY_DEVICE_INLINE double
total_energy(const EOS &eos, double rho, double v1, double v2, double p);

/// Are the values obviously unphysical?
/** The variable `u` is considered implausible if any component is not even
 *  a real number or if the density, or energy is non-positive.
 *
 *  @param u variable to check
 *  @return true if clearly unphysical.
 */
ANY_DEVICE_INLINE bool notplausible(const euler_var_t &u);

/// Could these values be physical (at first glace)?
/** The variable `u` is considered plausible if all its components are real
 *  numbers and furthermore, density and energy are (strictly) positive.
 *
 *  @param u variable to check
 *  @return false if clearly unphysical.
 */
ANY_DEVICE_INLINE bool isplausible(const euler_var_t &u);

template <class Model>
class euler_like {
public:
  static constexpr bool value = false;
};

template <>
class euler_like<Euler> {
public:
  static constexpr bool value = true;
};

} // namespace zisa

#endif /* end of include guard: EULER_H_OTPMLDSW */
