#include <catch/catch.hpp>
#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/weno_ao.hpp>

TEST_CASE("select stencil") {

  auto grid = zisa::load_gmsh("grids/convergence/unit_square_2.msh");
  SECTION("central_stencil") {

    zisa::int_t i_cell = 42;
    zisa::int_t n_points = 12;

    const auto &x_cell = grid->cell_centers(i_cell);

    auto stencil = zisa::central_stencil(*grid, i_cell, n_points);
    REQUIRE(stencil.size() == n_points);

    auto dist = [&x_cell, &grid](zisa::int_t i) {
      auto x = grid->cell_centers(i);
      return zisa::norm(x - x_cell);
    };

    auto is_in_stencil = [&stencil](zisa::int_t i) {
      return std::find(stencil.begin(), stencil.end(), i) != stencil.end();
    };

    std::sort(stencil.begin(),
              stencil.end(),
              [&dist](zisa::int_t i, zisa::int_t j) { return dist(i) < dist(j); });

    auto r_max = dist(stencil.back());

    for (zisa::int_t i = 0; i < grid->cell_centers.shape(0); ++i) {
      REQUIRE(((dist(i) > r_max) || is_in_stencil(i)));
    }
  }
}
