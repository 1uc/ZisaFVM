#include <memory>

#include <zisa/boundary/no_boundary_condition.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_integration_factory.hpp>
#include <zisa/ode/time_keeper_factory.hpp>

namespace zisa {

NumericalExperiment::NumericalExperiment(const InputParameters &params)
    : params(params),
      grid(choose_grid()),
      file_name_generator(choose_file_name_generator()) {}

void NumericalExperiment::run() { do_run(); }
void NumericalExperiment::post_process() { do_post_process(); }

std::shared_ptr<Grid> NumericalExperiment::choose_grid() {
  auto grid = load_gmsh(params["grid"]["file"]);
  return grid;
}

std::shared_ptr<FileNameGenerator>
NumericalExperiment::choose_file_name_generator() {
  return make_file_name_generator(params["io"]["filename"]);
}

void NumericalExperiment::write_grid() const {

// extra scope so writer goes out of scope and closes.
auto writer = HDF5SerialWriter("grid.h5");
save(writer, *grid);
}

void NumericalExperiment::do_run() {
  write_grid();

  std::cout << " --- Grid ---------- \n";
  std::cout << grid->str() << "\n";

  auto u0 = choose_initial_conditions();
  auto time_loop = choose_time_loop();

  auto u1 = (*time_loop)(u0);

  do_post_run(u1);
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
  auto instantaneous_physics = choose_instantaneous_physics();
  auto sanity_check = choose_sanity_check();
  auto visualization = choose_visualization();
  auto cfl_condition = choose_cfl_condition();

  return std::make_shared<TimeLoop>(time_integration,
                                    instantaneous_physics,
                                    simulation_clock,
                                    cfl_condition,
                                    sanity_check,
                                    visualization);
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

std::shared_ptr<RateOfChange> NumericalExperiment::aggregate_rates_of_change(
    const std::vector<std::shared_ptr<RateOfChange>>
        &physical_rates_of_change) {

  auto zero_change = std::make_shared<ZeroRateOfChange>();
  auto flux_bc = choose_flux_bc();

  auto roc = std::make_shared<SumRatesOfChange>();
  roc->add_term(zero_change);

  for (const auto &phys_roc : physical_rates_of_change) {
    roc->add_term(phys_roc);
  }

  roc->add_term(flux_bc);

  return roc;
}

std::shared_ptr<InstantaneousPhysics>
NumericalExperiment::choose_instantaneous_physics() {
  return std::make_shared<NoInstantaneousPhysics>();
}

} // namespace zisa
