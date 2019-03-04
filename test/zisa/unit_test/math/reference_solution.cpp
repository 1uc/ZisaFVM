#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/reference_solution.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>

namespace zisa {
template <class euler_t>
std::shared_ptr<AllVariables> set_initial_conditions(const Grid &grid,
                                                     const euler_t &euler) {

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
                               + euler.eos.internal_energy(RhoP{rho, p})};
  }

  return all_vars;
}

}

TEST_CASE("ReferenceSolution; basic API", "[math]") {
  auto fine_grid = zisa::load_gmsh("grids/convergence/unit_square_4.msh");
  auto coarse_grid = zisa::load_gmsh("grids/convergence/unit_square_0.msh");

  auto euler = zisa::make_default_euler();
  auto all_vars_ref = set_initial_conditions(*fine_grid, euler);

  auto ref = zisa::EulerReferenceSolution<
      zisa::IsentropicEquilibrium<decltype(euler)::eos_t,
                                  decltype(euler)::gravity_t>>(
      fine_grid, all_vars_ref, {euler.eos, euler.gravity, 4});
}
