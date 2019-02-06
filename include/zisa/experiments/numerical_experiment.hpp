#ifndef NUMERICAL_EXPERIMENT_H_349CI
#define NUMERICAL_EXPERIMENT_H_349CI

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/core/time_loop.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/triangular_rule.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/sanity_check.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_integration.hpp>

namespace zisa {

class NumericalExperiment {
public:
  NumericalExperiment(const InputParameters &params);
  virtual ~NumericalExperiment() = default;

  void run();

protected:
  virtual void do_run();

  virtual std::shared_ptr<Grid> choose_grid();

  /// Total rate of change in the system.
  /** See also:
   *    `aggregate_rates_of_change`
   **/
  virtual std::shared_ptr<RateOfChange> choose_rate_of_change() = 0;

  /// Combine the physical fluxes into the total rate of change.
  /** This will initialize a `SumRatesOfChange` as follows:
   *    - `ZeroRateOfChange` so that everyone can increment.
   *    - all physical rates of change.
   *    - flux boundary conditions.
   **/
  virtual std::shared_ptr<RateOfChange>
  aggregate_rates_of_change(const std::vector<std::shared_ptr<RateOfChange>>
                                &physical_rates_of_change);

  virtual std::shared_ptr<SanityCheck> choose_sanity_check() = 0;
  virtual std::shared_ptr<Visualization> choose_visualization() = 0;
  virtual std::shared_ptr<AllVariables> choose_initial_conditions() = 0;
  virtual std::shared_ptr<CFLCondition> choose_cfl_condition() = 0;
  virtual AllVariablesDimensions choose_all_variable_dims() = 0;
  virtual std::shared_ptr<RateOfChange> choose_flux_bc() = 0;

  virtual EdgeRule choose_edge_rule();
  virtual TriangularRule choose_volume_rule();

  virtual std::shared_ptr<TimeIntegration> choose_time_integration();
  virtual std::shared_ptr<SimulationClock> choose_simulation_clock();
  virtual std::shared_ptr<BoundaryCondition> choose_boundary_condition();
  virtual std::shared_ptr<TimeLoop> choose_time_loop();

protected:
  InputParameters params;

  std::shared_ptr<Grid> grid;
};

} // namespace zisa

#endif
