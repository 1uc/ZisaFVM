#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/experiments/numerical_experiment_factory.hpp>
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

  auto ic_ = PolytropeIC(euler);
  auto ic = [&eos, &ic_, amp, width](const auto &x) {
    // Avoid any silent conversion when the return type changes.
    RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;
    double p = p_eq * (1 + amp * exp(-zisa::pow<2>(r / width)));

    return eos.cvars(RhoP{rho_eq, p});
  };

  auto &u0 = all_variables->cvars;
  zisa::for_each(triangles(*grid),
                 [&qr, &u0, &ic](int_t i, const Triangle &tri) {
                   u0(i) = average(qr, ic, tri);
                 });

  return all_variables;
}

std::shared_ptr<AllVariables> JankaBump::compute_initial_conditions() {
  double amp = params["experiment"]["initial_conditions"]["amplitude"];
  double width = params["experiment"]["initial_conditions"]["width"];

  auto all_variables = compute_initial_conditions(amp, width);
  auto steady_state = compute_initial_conditions(0.0, width);

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
