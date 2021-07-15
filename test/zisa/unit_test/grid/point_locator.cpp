// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <zisa/grid/point_locator.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

TEST_CASE("PointLocator; basic API", "[grid][point_locator]") {
  auto grid_names = std::vector<std::string>{};

  for (int i = 0; i < 4; ++i) {
    grid_names.push_back(zisa::TestGridFactory::unit_square(i));
  }

  for (const auto &grid_name : grid_names) {
    auto grid = zisa::load_grid(grid_name);
    auto locator = make_point_locator(grid);

    auto n_cells = grid->n_cells;

    for (zisa::int_t i = 0; i < n_cells; ++i) {
      auto i_cell = locator->locate(grid->cell_centers(i));

      REQUIRE(i_cell == i);
    }
  }
}
