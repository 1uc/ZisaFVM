// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/loops/for_each.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/model/grid_variables.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

namespace zisa {

TEST_CASE("for_each; saxpy-like update") {
  auto grid = zisa::load_grid(zisa::TestGridFactory::unit_square(0));

  zisa::int_t n_vars = 3;
  zisa::int_t n_cells = grid->n_cells;

  auto u0 = zisa::GridVariables(zisa::shape_t<2>{n_cells, n_vars});
  auto u1 = zisa::GridVariables(zisa::shape_t<2>{n_cells, n_vars});

  zisa::fill(u1, -1.0);

  auto f = [](zisa::int_t i, zisa::int_t k) { return double(100 * i + k); };

  for (zisa::int_t i = 0; i < n_cells; ++i) {
    for (zisa::int_t k = 0; k < n_vars; ++k) {
      u0(i, k) = f(i, k);
    }
  }

  auto ff = [f, &u1, &u0](zisa::int_t i, const zisa::Cell &) {
    for (zisa::int_t k = 0; k < u0.shape(1); ++k) {
      u1(i, k) = f(i, k) + u0(i, k);
    }
  };

  zisa::for_each(zisa::cells(*grid), ff);

  for (zisa::int_t i = 0; i < n_cells; ++i) {
    for (zisa::int_t k = 0; k < n_vars; ++k) {
      auto approx = u1(i, k);
      auto exact = 2.0 * f(i, k);

      INFO(string_format(
          "(%d, %d): %e - %e = %e", i, k, approx, exact, approx - exact));
      REQUIRE(zisa::almost_equal(approx, exact, 1e-10));
    }
  }
}

}
