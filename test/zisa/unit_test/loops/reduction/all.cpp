#include <zisa/loops/reduction/all.hpp>

#include <zisa/model/grid_variables.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("all; basic examples") {
  zisa::int_t n_vars = 3;
  zisa::int_t n_cells = 10;

  auto grid_vars = zisa::GridVariables(zisa::shape_t<2>{n_cells, n_vars});

  auto is_valid
      = [](zisa::int_t, const auto &u) { return u[0] * u[1] * u[2] > 0; };

  SECTION("all true") {
    for (zisa::int_t i = 0; i < n_cells; ++i) {
      grid_vars(i, 0) = i + 2;
      grid_vars(i, 1) = i + 3;
      grid_vars(i, 2) = i + 2;
    }

    bool expected = true;
    bool observed = zisa::reduce::all(cell_const_range(grid_vars), is_valid);
    REQUIRE(observed == expected);
  }

  SECTION("some false") {
    for (zisa::int_t i = 0; i < n_cells; ++i) {
      grid_vars(i, 0) = i + 2;
      grid_vars(i, 1) = -double(i);
      grid_vars(i, 2) = i + 2;
    }

    bool expected = false;
    bool observed = zisa::reduce::all(cell_const_range(grid_vars), is_valid);
    REQUIRE(observed == expected);
  }
}