#include <zisa/boundary/frozen_boundary_condition.hpp>
#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/experiments/polytrope.hpp>
#include <zisa/parallelization/omp.h>

namespace zisa {

std::shared_ptr<AllVariables> Polytrope::compute_initial_conditions() {
  double amp = params["experiment"]["initial_conditions"]["amplitude"];
  double width = params["experiment"]["initial_conditions"]["width"];

  auto all_variables = compute_initial_conditions(amp, width);
  auto steady_state = compute_initial_conditions(0.0, width);

  auto vis = choose_visualization();
  vis->steady_state(*steady_state);

  return all_variables;
}

std::shared_ptr<AllVariables>
Polytrope::compute_initial_conditions(double amp, double width) {
  auto dims = choose_all_variable_dims();
  auto all_variables = std::make_shared<AllVariables>(dims);
  auto grid = choose_grid();

  auto qr = choose_volume_rule();
  const auto &eos = euler->eos;

  // FIXME revert
  auto ic_ = GeneralPolytropeIC(euler, {1.0, 1.0});
  auto ic = [&eos, &ic_, amp, width](const auto &x) {
    // Avoid any silent conversion when the return type changes.
    RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;
    double p = p_eq * (1 + amp * exp(-zisa::pow<2>(r / width)));

    return eos.cvars(RhoP{rho_eq, p});
  };

  auto &u0 = all_variables->cvars;
  zisa::for_each(cells(*grid), [&u0, &ic](int_t i, const Cell &cell) {
    u0(i) = average(cell, ic);
  });

  return all_variables;
}

std::function<bool(const Grid &grid, int_t i)>
Polytrope::boundary_mask() const {

  return [](const Grid &grid, int_t i) {
    // FIXME hard-coded value.
    double r_crit = 0.5;
    return zisa::norm(grid.cell_centers[i]) > r_crit;
  };
}

std::function<std::shared_ptr<Grid>(const std::string &, int_t)>
Polytrope::choose_grid_factory() {
  return [this](const std::string &filename, int_t quad_deg) {
    auto grid = load_grid(filename, quad_deg);
    enforce_cell_flags(*grid);
    return grid;
  };
}

std::shared_ptr<AllVariables> JankaBump::compute_initial_conditions() {
  double amp = params["experiment"]["initial_conditions"]["amplitude"];
  double width = params["experiment"]["initial_conditions"]["width"];

  auto all_variables = compute_initial_conditions(amp, width);
  auto steady_state = compute_initial_conditions(0.0, width);

  LOG_ERR("This needs to be updated.")
  auto fng = choose_file_name_generator();
  auto writer = HDF5SerialWriter(fng->steady_state_filename);
  save(writer, *steady_state, all_labels<euler_var_t>());

  return all_variables;
}

std::shared_ptr<AllVariables>
JankaBump::compute_initial_conditions(double amp, double width) {
  auto dims = choose_all_variable_dims();
  auto all_variables = std::make_shared<AllVariables>(dims);
  auto grid = choose_grid();

  auto qr = choose_volume_rule();
  const auto &eos = euler->eos;

  double rho_center = params["euler"]["gravity"]["rhoC"];
  double K_center = params["euler"]["gravity"]["K"];
  auto rhoK_center = RhoEntropy{rho_center, K_center};

  auto ic_ = GeneralPolytropeIC(euler, rhoK_center);
  auto ic = [&eos, &ic_, amp, width](const auto &x) {
    // Avoid any silent conversion when the return type changes.
    RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;
    double p = p_eq * (1 + amp * exp(-zisa::pow<2>(r / width)));

    return eos.cvars(RhoP{rho_eq, p});
  };

  auto &u0 = all_variables->cvars;

  auto n_cells = grid->n_cells;
#if ZISA_HAS_OPENMP == 1
#pragma omp parallel for ZISA_OMP_FOR_SCHEDULE_DEFAULT
#endif
  for (int_t i = 0; i < n_cells; ++i) {
    u0(i) = average(grid->cells(i), ic);
  }

  return all_variables;
}

} // namespace zisa
