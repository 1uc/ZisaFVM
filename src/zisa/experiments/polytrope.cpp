#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/experiments/polytrope.hpp>
#include <zisa/parallelization/omp.h>

namespace zisa {

std::shared_ptr<AllVariables> Polytrope::choose_initial_conditions() {
  double amp = params["experiment"]["initial_conditions"]["amplitude"];
  double width = params["experiment"]["initial_conditions"]["width"];

  auto all_variables = choose_initial_conditions(amp, width);
  auto steady_state = choose_initial_conditions(0.0, width);

  auto writer = HDF5SerialWriter(file_name_generator->steady_state_filename);
  save(writer, *steady_state, all_labels<euler_var_t>());

  return all_variables;
}

std::shared_ptr<AllVariables>
Polytrope::choose_initial_conditions(double amp, double width) {
  auto dims = choose_all_variable_dims();
  auto all_variables = std::make_shared<AllVariables>(dims);

  auto qr = choose_volume_rule();
  const auto &eos = euler.eos;

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

  auto n_cells = grid->n_cells;
#pragma omp parallel for ZISA_OMP_FOR_SCHEDULE_DEFAULT
  for (int_t i = 0; i < n_cells; ++i) {
    auto tri = grid->triangle(i);
    u0(i) = average(qr, ic, tri);
  }

  return all_variables;
}

} // namespace zisa
