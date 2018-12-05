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
     const = 0;

  /// Short self-documenting string.
  virtual std::string str() const = 0;
};

/// Convenient way to sum individual terms of the rate of change.
class SumRatesOfChange : public RateOfChange {
public:
  SumRatesOfChange() = default;
  SumRatesOfChange(std::vector<std::shared_ptr<RateOfChange>> init_list);

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double t) const override;

  void add_term(const std::shared_ptr<RateOfChange> &rate);

  void remove_all_terms();

  /// Short self-documenting string.
  virtual std::string str() const override;

private:
  std::vector<std::shared_ptr<RateOfChange>> rates_of_change;
};

/// Set the tendency buffer to zero.
class ZeroRateOfChange : public RateOfChange {
public:
  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double t) const override;

  virtual std::string str() const override;
};

} // namespace zisa
#endif /* end of include guard: RATE_OF_CHANGE_H_ES9JNRXD */
