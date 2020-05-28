#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/reference_solution.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

namespace zisa {
template <class euler_t>
std::shared_ptr<AllVariables> initial_conditions(const Grid &grid,
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
  auto fine_grid = zisa::load_grid(zisa::TestGridFactory::unit_square(3), 4);
  auto coarse_gridnames
      = std::vector<std::string>{zisa::TestGridFactory::unit_square(0),
                                 zisa::TestGridFactory::unit_square(1)};

  auto euler = zisa::make_default_euler();
  auto all_vars_ref = initial_conditions(*fine_grid, euler);

  using eq_t = zisa::NoEquilibrium;
  auto eq = eq_t{};

  using scaling_t = zisa::UnityScaling;
  auto scaling = scaling_t{};

  auto ref = zisa::EulerReferenceSolution<eq_t, scaling_t>(
      fine_grid, all_vars_ref, eq, scaling);

  zisa::int_t k_var = 0;

  for (const auto &coarse_gridname : coarse_gridnames) {
    auto coarse_grid = zisa::load_grid(coarse_gridname, 4);
    auto n_cells = coarse_grid->n_cells;

    auto approx = ref.average(*coarse_grid);
    auto exact = initial_conditions(*coarse_grid, euler);

    double l1_err = 0.0;
    for (zisa::int_t i = 0; i < n_cells; ++i) {

      double q_approx = approx->cvars(i, k_var);
      double q_exact = exact->cvars(i, k_var);

      l1_err += coarse_grid->volumes(i) * zisa::abs(q_approx - q_exact);

      double r_max = zisa::circum_radius(coarse_grid->triangle(i));

      double dq_max = r_max * 1.0;
      double q_min = q_exact - dq_max;
      double q_max = q_exact + dq_max;

      INFO(string_format("[%d] q = %.3e, dq = %.3e", i, q_exact, dq_max));
      REQUIRE(q_approx <= q_max);

      INFO(string_format("[%d] q = %.3e, dq = %.3e", i, q_exact, dq_max));
      REQUIRE(q_approx >= q_min);
    }
  }
}
