// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/reconstruction/lsq_solver.hpp>

#include <zisa/reconstruction/assemble_weno_ao_matrix.hpp>
#include <zisa/reconstruction/lsq_solver_family.hpp>
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

TEST_CASE("LSQSolver; assemble_weno_ao_matrix", "[lsq][3d]") {
  auto params = zisa::StencilFamilyParams(
      {3, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5, 1.5});

  auto quad_deg = max_order(params);
  auto grid = zisa::load_grid(zisa::TestGridFactory::unit_cube(0), quad_deg);

  auto n_cells = grid->n_cells;
  for (zisa::int_t i_cell = 0; i_cell < n_cells; ++i_cell) {
    auto stencils = zisa::StencilFamily(*grid, i_cell, params);

    double magic_value = -12345.6;

    for (const auto &stencil : stencils) {
      if (stencil.order() == 1) {
        continue;
      }

      auto A = zisa::allocate_weno_ao_matrix(*grid, stencil);
      A = Eigen::MatrixXd::Constant(A.rows(), A.cols(), magic_value);

      auto A_ref = Eigen::Ref<Eigen::MatrixXd>(A);

      zisa::assemble_weno_ao_matrix(A_ref, *grid, stencil);

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
