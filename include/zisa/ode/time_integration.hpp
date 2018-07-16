/* Explicit time stepper interface.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-02-01
 */
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
   *  @param current_state
   *      Values of the variables at time `t`.
   *  @param t  current time
   *  @param dt  current time step
   */
  virtual std::shared_ptr<AllVariables>
  compute_step(const std::shared_ptr<AllVariables> &u0, double t, double dt)
      = 0;

  /// Apply the boundary conditions of the PDE.
  virtual void boundary_condition(AllVariables &u0, double t) const = 0;

  /// Pick a stable time step smaller than `dt`.
  /** @param all_variables  current state prognostic variables.
   */
  virtual double pick_time_step(const AllVariables &all_variables) const = 0;

  /// Pick a stable time step smaller than `dt`.
  /** @param all_variables  current state prognostic variables.
   *  @param dt  proposed time-step.
   */
  virtual double pick_time_step(const AllVariables &all_variables,
                                double dt) const = 0;

  /// Self-documenting string.
  virtual std::string str(int indent = 0) const = 0;
};

} // namespace zisa
#endif /* end of include guard: TIME_INTEGRATION_H_N9QBMXXY */