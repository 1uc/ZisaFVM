// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef NUMERICAL_EXPERIMENT_H_349CI
#define NUMERICAL_EXPERIMENT_H_349CI

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/cli/input_parameters.hpp>
#include <zisa/fvm_loops/time_loop.hpp>
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
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/reconstruction/stencil_family_params.hpp>

namespace zisa {

class NumericalExperiment {
public:
  virtual ~NumericalExperiment() = default;

  void run();
  void post_process();

protected:
  virtual void do_run() = 0;
  virtual void do_post_process() = 0;
};

class InvalidNumericalExperiment : public NumericalExperiment {
public:
  InvalidNumericalExperiment(std::string error_message)
      : error_message(std::move(error_message)) {}

protected:
  virtual void do_run() override { LOG_ERR(error_message); }
  virtual void do_post_process() override { LOG_ERR(error_message); }

private:
  std::string error_message;
};

class TypicalNumericalExperiment : public NumericalExperiment {
public:
  explicit TypicalNumericalExperiment(const InputParameters &params);

protected:
  virtual void write_grid();
  virtual void do_run() override;
  virtual void do_post_run(const std::shared_ptr<AllVariables> &u1) = 0;

  bool is_restart() const;

  std::shared_ptr<Grid> choose_grid() const;
  virtual std::shared_ptr<Grid> compute_grid() const;

  virtual std::shared_ptr<array<StencilFamily, 1>> choose_stencils() const;

  virtual std::shared_ptr<array<StencilFamily, 1>>
  compute_stencils(const Grid &grid) const;

  virtual StencilFamilyParams choose_stencil_params() const;

  virtual void print_grid_info();

  std::shared_ptr<FileNameGenerator> choose_file_name_generator();
  virtual std::shared_ptr<FileNameGenerator> compute_file_name_generator();

  std::pair<std::shared_ptr<AllVariables>, std::shared_ptr<AllVariables>>
  choose_initial_conditions();

  virtual std::pair<std::shared_ptr<AllVariables>,
                    std::shared_ptr<AllVariables>>
  compute_initial_conditions() = 0;
  virtual std::pair<std::shared_ptr<AllVariables>,
                    std::shared_ptr<AllVariables>>
  load_initial_conditions() = 0;

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

  virtual std::shared_ptr<Visualization> choose_visualization();
  virtual std::shared_ptr<Visualization> compute_visualization() = 0;

  virtual std::shared_ptr<CFLCondition> choose_cfl_condition() = 0;
  virtual AllVariablesDimensions choose_all_variable_dims() = 0;
  virtual int_t choose_n_avars();
  virtual std::shared_ptr<RateOfChange> choose_flux_bc() = 0;

  virtual EdgeRule choose_edge_rule();
  virtual TriangularRule choose_volume_rule();

  int_t choose_volume_deg() const;
  int_t choose_edge_deg() const;
  int_t choose_moments_deg() const;
  QRDegrees choose_qr_degrees() const;

  virtual std::shared_ptr<TimeIntegration> choose_time_integration();

  virtual std::shared_ptr<SimulationClock> choose_simulation_clock();
  virtual std::shared_ptr<SimulationClock> compute_simulation_clock();

  std::shared_ptr<BoundaryCondition> choose_boundary_condition();
  virtual std::shared_ptr<BoundaryCondition> compute_boundary_condition();

  virtual std::shared_ptr<TimeLoop> choose_time_loop();
  virtual std::shared_ptr<ProgressBar> choose_progress_bar();

  virtual void enforce_cell_flags(Grid &grid) const;
  virtual std::function<bool(const Grid &, int_t)> boundary_mask() const;

  virtual void write_debug_output();
  virtual void write_stencils();

protected:
  InputParameters params;
  std::shared_ptr<FileNameGenerator> file_name_generator_;

protected:
  mutable std::shared_ptr<Grid> grid_ = nullptr;
  mutable std::shared_ptr<Visualization> visualization_ = nullptr;
  mutable std::shared_ptr<BoundaryCondition> boundary_condition_ = nullptr;
  mutable std::shared_ptr<SimulationClock> simulation_clock_ = nullptr;
  mutable std::shared_ptr<array<StencilFamily, 1>> stencils_ = nullptr;
  mutable std::shared_ptr<AllVariables> all_vars_ = nullptr;
  mutable std::shared_ptr<AllVariables> steady_state_ = nullptr;

  time_stamp_t t_start_ = current_time_stamp();
};

} // namespace zisa

#endif
