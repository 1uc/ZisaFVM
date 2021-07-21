// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/grid/cell_range.hpp>

#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

TEST_CASE("CellRange; API", "[ranges][runme]") {
  auto grid = zisa::load_grid(zisa::TestGridFactory::unit_square(0));

  for (auto [i, cell] : zisa::cells(*grid)) {
    REQUIRE(grid->cells(i) == cell);
  }
}