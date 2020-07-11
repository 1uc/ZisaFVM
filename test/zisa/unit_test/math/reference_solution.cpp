#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/reference_solution.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

namespace zisa {
template <class EOS>
std::shared_ptr<AllVariables> initial_conditions(const Grid &grid,
                                                 const EOS &eos) {

  int_t n_cells = grid.n_cells;
  int_t n_cvars = 5;
  int_t n_avars = 0;

  auto all_vars = std::make_shared<zisa::AllVariables>(
      zisa::AllVariablesDimensions{n_cells, n_cvars, n_avars});

  auto &cvars = all_vars->cvars;

  for (int_t i = 0; i < n_cells; ++i) {
    const auto &x = grid.cell_centers(i);

    double rho = 1.0 + 0.5 * (zisa::sin(x[0] + x[1]) + zisa::cos(x[1] * x[0]));
    double vx = 1.0 + 0.5 * (zisa::sin(x[0] * x[1]) + zisa::cos(x[1] * x[0]));
    double vy = 1.0 + 0.5 * (zisa::sin(x[1] * x[1]) * zisa::cos(x[1] * x[0]));
    double p = 20.0;

    cvars(i) = euler_var_t{rho,
                           rho * vx,
                           rho * vy,
                           0.0,
                           0.5 * rho * (vx * vx + vy * vy)
                               + eos.internal_energy(RhoP{rho, p})};
  }

  return all_vars;
}

}
