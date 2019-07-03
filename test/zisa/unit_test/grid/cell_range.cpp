#include <zisa/grid/cell_range.hpp>

#include <zisa/testing/testing_framework.hpp>

TEST_CASE("CellRange; API", "[ranges][runme]") {
  auto grid = zisa::load_gmsh("grids/convergence/unit_square_0.msh", 1);

  for (auto [i, cell] : zisa::cells(*grid)) {
    REQUIRE(grid->cells(i) == cell);
  }
}