#include <memory>

#include <zisa/boundary/no_boundary_condition.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_integration_factory.hpp>
#include <zisa/ode/time_keeper_factory.hpp>

namespace zisa {

NumericalExperiment::NumericalExperiment(const InputParameters &params)
    : params(params) {}

void NumericalExperiment::run() { do_run(); }

std::shared_ptr<Grid> NumericalExperiment::choose_grid() {
  return load_gmsh(params["grid"]["file"]);
}

void NumericalExperiment::do_run() {
  grid = choose_grid();

  auto u0 = choose_initial_conditions();
  auto time_loop = choose_time_loop();

  (*time_loop)(u0);
}

std::shared_ptr<SimulationClock>
NumericalExperiment::choose_simulation_clock() {

  auto time_keeper_params = TimeKeeperParameters(params);
  auto plotting_params = PlottingStepsParameters(params);

  auto time_keeper = make_time_keeper(time_keeper_params);
  auto plotting_steps
      = make_plotting_steps(plotting_params, time_keeper_params);

  return std::make_shared<SerialSimulationClock>(time_keeper, plotting_steps);
}

std::shared_ptr<TimeLoop> NumericalExperiment::choose_time_loop() {
  auto simulation_clock = choose_simulation_clock();
  auto time_integration = choose_time_integration();
  auto sanity_check = choose_sanity_check();
  auto visualization = choose_visualization();
  auto cfl_condition = choose_cfl_condition();

  return std::make_shared<TimeLoop>(
      time_integration, simulation_clock, cfl_condition, sanity_check, visualization);
}

std::shared_ptr<BoundaryCondition>
NumericalExperiment::choose_boundary_condition() {
  return std::make_shared<NoBoundaryCondition>();
}

EdgeRule NumericalExperiment::choose_edge_rule() {
  LOG_ERR_IF(!has_key(params, "quadrature"), "Missing section 'quadrature'.");
  LOG_ERR_IF(!has_key(params["quadrature"], "edge"), "Missing element 'edge'.");

  return cached_edge_quadrature_rule(params["quadrature"]["edge"]);
}

TriangularRule NumericalExperiment::choose_volume_rule() {
  LOG_ERR_IF(!has_key(params, "quadrature"), "Missing section 'quadrature'.");
  LOG_ERR_IF(!has_key(params["quadrature"], "volume"),
             "Missing element 'volume'.");

  return cached_triangular_quadrature_rule(params["quadrature"]["volume"]);
}

std::shared_ptr<GlobalReconstruction<CWENO_AO>>
NumericalExperiment::choose_reconstruction() {
  LOG_ERR_IF(!has_key(params, "reconstruction"),
             "Missing section 'reconstruction'.");
  auto rc_params = params["reconstruction"];

  auto hybrid_weno_params
      = HybridWENO_Params(StencilFamilyParams(rc_params["orders"],
                                              rc_params["biases"],
                                              rc_params["overfit_factors"]),
                          rc_params["linear_weights"],
                          rc_params["smoothness_indicator"]["epsilon"],
                          rc_params["smoothness_indicator"]["exponent"]);

  auto dims = choose_all_variable_dims();

  return std::make_shared<GlobalReconstruction<CWENO_AO>>(
      grid, hybrid_weno_params, dims.n_cvars);
}

std::shared_ptr<TimeIntegration>
NumericalExperiment::choose_time_integration() {
  LOG_ERR_IF(!has_key(params, "ode"),
             "Missing section 'ode' in the config file.");

  LOG_ERR_IF(!has_key(params["ode"], "solver"),
             "Missing key 'solver' in the config file.");

  std::string desc = params["ode"]["solver"];

  auto bc = choose_boundary_condition();
  auto dims = choose_all_variable_dims();
  auto rate_of_change = choose_rate_of_change();

  return make_time_integration(desc, rate_of_change, bc, dims);
}

} // namespace zisa
