#ifndef NUMERICAL_EXPERIMENT_H_349CI
#define NUMERICAL_EXPERIMENT_H_349CI

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/core/time_loop.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/io/file_name_generator.hpp>
#include <zisa/io/visualization.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/triangular_rule.hpp>
#include <zisa/model/cfl_condition.hpp>
#include <zisa/model/instantaneous_physics.hpp>
#include <zisa/model/sanity_check.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/step_rejection.hpp>
#include <zisa/ode/time_integration.hpp>

namespace zisa {

class NumericalExperiment {
public:
  explicit NumericalExperiment(const InputParameters &params);
  virtual ~NumericalExperiment() = default;

  void run();
  void post_process();

protected:
  virtual void write_grid() const;
  virtual void do_run();
  virtual void do_post_run(const std::shared_ptr<AllVariables> &u1) = 0;
  virtual void do_post_process() = 0;

  bool is_restart() const { return has_key(params, "restart"); }

  std::shared_ptr<Grid> choose_grid();
  std::shared_ptr<FileNameGenerator> choose_file_name_generator();

  std::shared_ptr<AllVariables> choose_initial_conditions();
  virtual std::shared_ptr<AllVariables> compute_initial_conditions() = 0;
  virtual std::shared_ptr<AllVariables> load_initial_conditions() = 0;

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

  virtual std::shared_ptr<InstantaneousPhysics> choose_instantaneous_physics();
  virtual std::shared_ptr<StepRejection> choose_step_rejection();
  virtual std::shared_ptr<SanityCheck> choose_sanity_check() = 0;
  virtual std::shared_ptr<Visualization> choose_visualization() = 0;
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
  std::shared_ptr<FileNameGenerator> file_name_generator;
};

} // namespace zisa

#endif
