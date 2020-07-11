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
  auto local_eos = choose_local_eos();
  const auto &eos = (*local_eos)(0);

  auto qr = choose_volume_rule();

  // FIXME revert
  auto ic_ = GeneralPolytropeIC(eos, gravity, {1.0, 1.0});
  auto ic = [&eos, &ic_, amp, width](const auto &x) {
    // Avoid any silent conversion when the return type changes.
    RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;
    double p = p_eq * (1 + amp * exp(-zisa::pow<2>(r / width)));

    return eos->cvars(RhoP{rho_eq, p});
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

} // namespace zisa
