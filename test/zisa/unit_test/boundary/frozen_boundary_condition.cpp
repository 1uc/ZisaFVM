// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <zisa/boundary/frozen_boundary_condition.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

TEST_CASE("FrozenBC; example", "[bc]") {
  auto grid_ = zisa::load_grid(zisa::TestGridFactory::unit_square(0));
  auto &grid = *grid_;

  auto n_cells = grid.n_cells;
  auto n_cvars = 5ul;
  auto n_avars = 0ul;

  auto I = std::vector<zisa::int_t>{0, 1, 4, 8, 16};
  for (auto i : I) {
    grid.cell_flags(i).interior = false;
    grid.cell_flags(i).ghost_cell = true;
  }

  auto dims = zisa::AllVariablesDimensions{n_cells, n_cvars, n_avars};
  auto all_vars = zisa::AllVariables(dims);
  auto &u = all_vars.cvars;

  zisa::fill(u, 1.0);

  for (auto i : I) {
    u(i, 0) = 100.0 * i + 0.0;
    u(i, 1) = 100.0 * i + 1.0;
    u(i, 2) = 100.0 * i + 2.0;
    u(i, 3) = 100.0 * i + 3.0;
    u(i, 4) = 100.0 * i + 4.0;
  }

  auto bc = std::make_shared<zisa::FrozenBC>(grid, all_vars);
  zisa::fill(u, 0.0);

  bc->apply(all_vars, /* t = */ 0.0);

  for (auto i : I) {
    REQUIRE(u(i, 0) == 100.0 * i + 0.0);
    REQUIRE(u(i, 1) == 100.0 * i + 1.0);
    REQUIRE(u(i, 2) == 100.0 * i + 2.0);
    REQUIRE(u(i, 3) == 100.0 * i + 3.0);
    REQUIRE(u(i, 4) == 100.0 * i + 4.0);
  }
}
