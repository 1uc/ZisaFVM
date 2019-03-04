#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/reference_solution.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_factory.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>

TEST_CASE("ReferenceSolution; basic API", "[math]") {
  auto fine_grid = zisa::load_gmsh("grid/convergence/unit_square_4.msh");
  auto coarse_grid = zisa::load_gmsh("grid/convergence/unit_square_0.msh");

  zisa::int_t n_cells = fine_grid->n_cells;
  zisa::int_t n_cvars = 5;
  zisa::int_t n_avars = 0;

  auto all_vars_ref = std::make_shared<zisa::AllVariables>(
      zisa::AllVariablesDimensions{n_cells, n_cvars, n_avars});

  auto euler = zisa::make_default_euler();

  auto ref = zisa::EulerReferenceSolution<
      zisa::IsentropicEquilibrium<decltype(euler)::eos_t,
                                  decltype(euler)::gravity_t>>(
      fine_grid, all_vars_ref, {euler.eos, euler.gravity, 4});
}
