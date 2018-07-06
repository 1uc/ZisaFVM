/* Interface for ODE rate of change.
 */
#ifndef RATE_OF_CHANGE_H_ES9JNRXD
#define RATE_OF_CHANGE_H_ES9JNRXD

#include "zisa/config.hpp"
#include "zisa/model/all_variables.hpp"

namespace zisa {
/// Abstract interface to compute the rate of change of an ODE.
/** A PDE can be discretized in two steps, first a spacial discretization,
 *  leaving an ODE with a very expensive RHS. The so called semi-discrete
 *  formulation.
 *
 *  This is the interface that computes this RHS and can be used in
 *  conjunction with `RungeKutta` or other time integrators.
 *
 *  @see TimeIntegration
 */
class RateOfChange {
public:
  virtual ~RateOfChange() = default;

  /// Compute rate of change.
  /** Compute the rate of change of the conserved variables.
   *
   *  @param[out] tendency
   *     The rate of change of the conserved variables.
   *  @param[in] current_state
   *     The current values of the conserved variables.
   *  @param[in] t
   *     Current time of the simulation.
   */
  virtual void
  compute(AllVariables &tendency, const AllVariables &current_state, double t)
      = 0;

  /// Pick a stable time-step.
  /** @param all_variables  current state prognostic variables.
   */
  virtual double pick_time_step(const AllVariables &all_variables) const = 0;

  /// Pick a stable time step smaller than `dt`.
  /** @param all_variables  current state prognostic variables.
   *  @param dt  proposed time-step.
   */
  virtual double pick_time_step(const AllVariables &all_variables,
                                double dt) const = 0;

  /// Short self-documenting string.
  virtual std::string str(int indent) const = 0;
};

/// Convenient way to sum individual terms of the rate of change.
class SumRatesOfChange : public RateOfChange {
public:
  SumRatesOfChange() = default;
  SumRatesOfChange(std::vector<std::shared_ptr<RateOfChange>> init_list);

  virtual ~SumRatesOfChange() = default;

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double t) override;

  void add_term(const std::shared_ptr<RateOfChange> &rate);

  void remove_all_terms();

  virtual double
  pick_time_step(const AllVariables &all_variables) const override;

  virtual double pick_time_step(const AllVariables &all_variables,
                                double dt) const override;

  /// Short self-documenting string.
  virtual std::string str(int indent) const override;

private:
  std::vector<std::shared_ptr<RateOfChange>> rates_of_change;
};

/// Set the tendency buffer to zero.
class ZeroRateOfChange : public RateOfChange {
public:
  ZeroRateOfChange(double dt_max);

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double t) override;

  virtual double
  pick_time_step(const AllVariables &all_variables) const override;

  virtual double pick_time_step(const AllVariables &all_variables,
                                double dt) const override;

  virtual std::string str(int indent) const override;

protected:
  double dt_max;
};

} // namespace zisa
#endif /* end of include guard: RATE_OF_CHANGE_H_ES9JNRXD */
