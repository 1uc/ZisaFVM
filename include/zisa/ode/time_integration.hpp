// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef TIME_INTEGRATION_H_N9QBMXXY
#define TIME_INTEGRATION_H_N9QBMXXY

#include "zisa/config.hpp"
#include "zisa/model/all_variables_fwd.hpp"

namespace zisa {

/// Interface of time integration schemes.
class TimeIntegration {
public:
  virtual ~TimeIntegration() = default;

  /// Evolve solution from t to t + dt.
  /** The caller must ensure that the time-step is sufficiently small
   *  to ensure that the integration is stable.
   *
   * Note: the update will not modify `u0`. Therefore, the returned smart
   *    pointer does not point to the same object as `u0`.
   *
   * Note: `RungeKutta` takes shared ownership of `u0` and therefore,
   *   can modify `u0` in a later call.
   *
   * These two comments imply that code such as
   *
   *     auto u1 = rk.compute_step(u0, t, dt);
   *     auto du = difference(*u0, *u1);
   *
   * is safe.
   *
   * WARNING: The two comments also imply that passing the same buffer
   *   twice is not safe.
   *
   *  @param current_state
   *      Values of the variables at time `t`.
   *  @param t  current time
   *  @param dt  current time step
   */
  virtual std::shared_ptr<AllVariables>
  compute_step(const std::shared_ptr<AllVariables> &u0, double t, double dt)
      = 0;

  /// Self-documenting string.
  virtual std::string str() const = 0;
};

} // namespace zisa
#endif /* end of include guard: TIME_INTEGRATION_H_N9QBMXXY */
