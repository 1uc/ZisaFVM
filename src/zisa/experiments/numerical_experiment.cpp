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
    : params(params), file_name_generator(choose_file_name_generator()) {}

void NumericalExperiment::run() { do_run(); }
void NumericalExperiment::post_process() { do_post_process(); }

std::shared_ptr<Grid> NumericalExperiment::choose_grid() const {
  if (grid_ == nullptr) {
    grid_ = compute_grid();
  }

  return grid_;
}

std::shared_ptr<Grid> NumericalExperiment::compute_grid() const {
  int_t quad_deg = params["quadrature"]["volume"];
  return load_grid(params["grid"]["file"], quad_deg);
}

std::shared_ptr<Grid> NumericalExperiment::choose_full_grid() const {
  if (full_grid_ == nullptr) {
    full_grid_ = compute_full_grid();
  }

  return full_grid_;
}

std::shared_ptr<Grid> NumericalExperiment::compute_full_grid() const {
  return choose_grid();
}

std::shared_ptr<FileNameGenerator>
NumericalExperiment::choose_file_name_generator() {
  auto fng = make_file_name_generator(params["io"]["filename"]);

  if (is_restart()) {
    fng->advance_to(std::string(params["restart"]["file"]));
  }

  return fng;
}

void NumericalExperiment::write_grid() {
  auto writer = HDF5SerialWriter("grid.h5");
  save(writer, *choose_grid());
}

void NumericalExperiment::do_run() {
  auto grid = choose_grid();

  if (!is_restart()) {
    write_grid();
  }

  print_grid_info();

  auto u0 = choose_initial_conditions();
  auto time_loop = choose_time_loop();

  auto u1 = (*time_loop)(u0);

  do_post_run(u1);
}

void NumericalExperiment::print_grid_info() {
  std::cout << " --- Grid ---------- \n";
  std::cout << choose_grid()->str() << "\n";
}

std::shared_ptr<AllVariables> NumericalExperiment::choose_initial_conditions() {
  if (is_restart()) {
    return load_initial_conditions();
  } else {
    return compute_initial_conditions();
  }
}

std::shared_ptr<SimulationClock>
NumericalExperiment::choose_simulation_clock() {

  auto time_keeper_params = TimeKeeperParameters(params);
  auto plotting_params = PlottingStepsParameters(params);

  auto time_keeper = make_time_keeper(time_keeper_params);
  auto plotting_steps
      = make_plotting_steps(plotting_params, time_keeper_params);

  auto simulation_clock
      = std::make_shared<SerialSimulationClock>(time_keeper, plotting_steps);

  if (is_restart()) {
    std::string datafile = params["restart"]["file"];

    auto reader = HDF5SerialReader(datafile);
    auto t = reader.read_scalar<double>("time");
    auto k = reader.read_scalar<int_t>("n_steps");

    simulation_clock->advance_to(t, k);
  }

  return simulation_clock;
}

std::shared_ptr<TimeLoop> NumericalExperiment::choose_time_loop() {
  auto simulation_clock = choose_simulation_clock();
  auto time_integration = choose_time_integration();
  auto instantaneous_physics = choose_instantaneous_physics();
  auto step_rejection = choose_step_rejection();
  auto sanity_check = choose_sanity_check();
  auto visualization = choose_visualization();
  auto cfl_condition = choose_cfl_condition();
  auto progress_bar = choose_progress_bar();

  return std::make_shared<TimeLoop>(time_integration,
                                    instantaneous_physics,
                                    step_rejection,
                                    simulation_clock,
                                    cfl_condition,
                                    sanity_check,
                                    visualization,
                                    progress_bar);
}

std::shared_ptr<BoundaryCondition>
NumericalExperiment::choose_boundary_condition() {
  return std::make_shared<NoBoundaryCondition>();
}

int_t NumericalExperiment::choose_volume_deg() const {
  return params["quadrature"]["volume"];
}

EdgeRule NumericalExperiment::choose_edge_rule() {
  return cached_edge_quadrature_rule(choose_edge_deg());
}

int_t NumericalExperiment::choose_edge_deg() const {
  return params["quadrature"]["edge"];
}

TriangularRule NumericalExperiment::choose_volume_rule() {
  return cached_triangular_quadrature_rule(choose_volume_deg());
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

std::shared_ptr<StepRejection> NumericalExperiment::choose_step_rejection() {
  return std::make_shared<RejectNothing>();
}

std::shared_ptr<ProgressBar> NumericalExperiment::choose_progress_bar() {
  return std::make_shared<SerialProgressBar>(1);
}

std::shared_ptr<Visualization> NumericalExperiment::choose_visualization() {
    if(visualization_ == nullptr) {
        visualization_ = compute_visualization();
    }

    return visualization_;
}

} // namespace zisa
