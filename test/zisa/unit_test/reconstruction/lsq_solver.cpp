#include <zisa/reconstruction/lsq_solver.hpp>

#include <zisa/reconstruction/lsq_solver_family.hpp>
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/testing/testing_framework.hpp>

namespace zisa {
Eigen::MatrixXd allocate_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

void assemble_weno_ao_matrix(Eigen::MatrixXd &A,
                             const Grid &grid,
                             const Stencil &stencil);
}

TEST_CASE("LSQSolver; assemble_weno_ao_matrix", "[lsq][3d]") {
  auto params = zisa::StencilFamilyParams(
      {3, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5, 1.5});

  auto grid = zisa::load_gmsh("grids/convergence/unit_cube_1.msh");

  auto n_cells = grid->n_cells;
  for (zisa::int_t i_cell = 0; i_cell < n_cells; ++i_cell) {
    auto stencils = zisa::StencilFamily(grid, i_cell, params);

    double magic_value = -12345.6;

    for (const auto &stencil : stencils) {
      if (stencil.order() == 1) {
        continue;
      }

      auto A = zisa::allocate_weno_ao_matrix(*grid, stencil);
      A = Eigen::MatrixXd::Constant(A.rows(), A.cols(), magic_value);

      zisa::assemble_weno_ao_matrix(A, *grid, stencil);

      for (Eigen::Index i = 0; i < A.rows(); ++i) {
        for (Eigen::Index j = 0; j < A.cols(); ++j) {
          INFO(string_format("[%d] A(%d, %d) = %f", i_cell, i, j, A(i, j)));
          REQUIRE(A(i, j) != magic_value);
          REQUIRE(zisa::isreal(A(i, j)));
        }
      }
    }
  }
}

