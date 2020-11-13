#include <memory>

#include <zisa/boundary/boundary_condition_factory.hpp>
#include <zisa/boundary/no_boundary_condition.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/memory/array_stencil_family.hpp>
#include <zisa/ode/simulation_clock.hpp>
#include <zisa/ode/time_integration_factory.hpp>
#include <zisa/ode/time_keeper_factory.hpp>

namespace zisa {

void NumericalExperiment::run() { do_run(); }
void NumericalExperiment::post_process() { do_post_process(); }

TypicalNumericalExperiment::TypicalNumericalExperiment(
    const InputParameters &params)
    : params(params) {}

std::shared_ptr<Grid> TypicalNumericalExperiment::choose_grid() const {
  if (grid_ == nullptr) {
    grid_ = compute_grid();
  }

  return grid_;
}

std::shared_ptr<Grid> TypicalNumericalExperiment::compute_grid() const {
  auto qr_degrees = choose_qr_degrees();
  auto grid = load_grid(params["grid"]["file"], qr_degrees);
  enforce_cell_flags(*grid);

  return grid;
}

std::shared_ptr<FileNameGenerator>
TypicalNumericalExperiment::choose_file_name_generator() {
  if (file_name_generator_ == nullptr) {
    file_name_generator_ = compute_file_name_generator();
  }

  return file_name_generator_;
}

std::shared_ptr<FileNameGenerator>
TypicalNumericalExperiment::compute_file_name_generator() {
  // needs to be implemented.
  assert(!is_restart());

  const auto &fn_params = params["io"]["filename"];

  auto fng = make_file_name_generator(
      "./", fn_params["stem"], fn_params["pattern"], fn_params["suffix"]);
  return fng;
}

void TypicalNumericalExperiment::write_grid() {
  auto writer = HDF5SerialWriter("grid.h5");
  save(writer, *choose_grid());
}

void TypicalNumericalExperiment::do_run() {
  auto grid = choose_grid();

  if (!is_restart()) {
    write_grid();
    write_debug_output();
  }

  print_grid_info();

  auto bc = choose_boundary_condition();
  auto u0 = choose_initial_conditions();
  bc->apply(*u0, /* t = */ 0.0);

  auto time_loop = choose_time_loop();

  auto u1 = (*time_loop)(u0);

  do_post_run(u1);
}

void TypicalNumericalExperiment::print_grid_info() {
  std::cout << " --- Grid ---------- \n";
  std::cout << choose_grid()->str() << "\n";
}

std::shared_ptr<AllVariables>
TypicalNumericalExperiment::choose_initial_conditions() {
  if (all_vars_ == nullptr) {
    if (is_restart()) {
      all_vars_ = load_initial_conditions();
    } else {
      all_vars_ = compute_initial_conditions();
    }
  }

  return all_vars_;
}

std::shared_ptr<SimulationClock>
TypicalNumericalExperiment::choose_simulation_clock() {

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

StencilFamilyParams TypicalNumericalExperiment::choose_stencil_params() const {
  const auto &rc_params = params["reconstruction"];
  return StencilFamilyParams(
      rc_params["orders"], rc_params["biases"], rc_params["overfit_factors"]);
}

std::shared_ptr<array<StencilFamily, 1>>
TypicalNumericalExperiment::choose_stencils() const {
  if (stencils_ == nullptr) {
    auto grid = choose_grid();
    stencils_ = compute_stencils(*grid);
  }

  return stencils_;
}

std::shared_ptr<array<StencilFamily, 1>>
TypicalNumericalExperiment::compute_stencils(const Grid &grid) const {
  assert(grid_ != nullptr);

  auto stencil_params = choose_stencil_params();

  return std::make_shared<array<StencilFamily, 1>>(
      compute_stencil_families(grid, stencil_params));
}

std::shared_ptr<TimeLoop> TypicalNumericalExperiment::choose_time_loop() {
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
TypicalNumericalExperiment::choose_boundary_condition() {
  if (boundary_condition_ == nullptr) {
    boundary_condition_ = compute_boundary_condition();
  }
  return boundary_condition_;
}

std::shared_ptr<BoundaryCondition>
TypicalNumericalExperiment::compute_boundary_condition() {
  auto grid = choose_grid();
  auto u0 = choose_initial_conditions();

  return make_boundary_condition(params, grid, u0);
}

int_t TypicalNumericalExperiment::choose_volume_deg() const {
  return params["quadrature"]["volume"];
}

EdgeRule TypicalNumericalExperiment::choose_edge_rule() {
  return cached_edge_quadrature_rule(choose_edge_deg());
}

int_t TypicalNumericalExperiment::choose_edge_deg() const {
  return params["quadrature"]["edge"];
}

int_t TypicalNumericalExperiment::choose_moments_deg() const {
  return params["quadrature"]["moments"];
}

QRDegrees TypicalNumericalExperiment::choose_qr_degrees() const {
  auto face_deg = choose_edge_deg();
  auto volume_deg = choose_volume_deg();
  auto moments_deg = choose_moments_deg();

  return QRDegrees{face_deg, volume_deg, moments_deg};
}

TriangularRule TypicalNumericalExperiment::choose_volume_rule() {
  return cached_triangular_quadrature_rule(choose_volume_deg());
}

std::shared_ptr<TimeIntegration>
TypicalNumericalExperiment::choose_time_integration() {
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

std::shared_ptr<RateOfChange>
TypicalNumericalExperiment::aggregate_rates_of_change(
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
TypicalNumericalExperiment::choose_instantaneous_physics() {
  return std::make_shared<NoInstantaneousPhysics>();
}

std::shared_ptr<StepRejection>
TypicalNumericalExperiment::choose_step_rejection() {
  return std::make_shared<RejectNothing>();
}

std::shared_ptr<ProgressBar> TypicalNumericalExperiment::choose_progress_bar() {
  return std::make_shared<SerialProgressBar>(1);
}

std::shared_ptr<Visualization>
TypicalNumericalExperiment::choose_visualization() {
  if (visualization_ == nullptr) {
    visualization_ = compute_visualization();
  }

  return visualization_;
}

void TypicalNumericalExperiment::write_debug_output() {
  if (has_key(params, "debug")) {
    if (params["debug"].value("stencils", false)) {
      write_stencils();
    }
  }
}

void TypicalNumericalExperiment::write_stencils() {
  const auto &stencils = choose_stencils();
  auto writer = HDF5SerialWriter("stencils.h5");
  save(writer, *stencils, "stencils");
}

void TypicalNumericalExperiment::enforce_cell_flags(Grid &grid) const {
  mask_ghost_cells(grid, boundary_mask());
}

std::function<bool(const Grid &, int_t)>
TypicalNumericalExperiment::boundary_mask() const {
  return [](const Grid &, int_t) { return false; };
}

int_t TypicalNumericalExperiment::choose_n_avars() { return 0; }

} // namespace zisa
