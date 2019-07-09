#include <numeric>
#include <zisa/testing/testing_framework.hpp>

#include <zisa/reconstruction/stencil.hpp>

TEST_CASE("select stencil", "[stencil][weno_ao]") {
  // FIXME remove argument
  auto grid = zisa::load_gmsh("grids/convergence/unit_square_2.msh", 0);
  zisa::int_t i_cell = 42;
  zisa::int_t n_points = 12;
  zisa::int_t n_cells = grid->n_cells;

  SECTION("central_stencil") {
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

    std::sort(
        stencil.begin(), stencil.end(), [&dist](zisa::int_t i, zisa::int_t j) {
          return dist(i) < dist(j);
        });

    auto r_max = dist(stencil.back());

    for (zisa::int_t i = 0; i < n_cells; ++i) {
      REQUIRE(((dist(i) > r_max) || is_in_stencil(i)));
    }
  }

  SECTION("biased_stencil; k = 0") {
    zisa::int_t k = 0;
    auto stencil = zisa::biased_stencil(*grid, i_cell, k, n_points);
    REQUIRE(stencil.size() == n_points);
  }

  SECTION("biased_stencil; k = 1") {
    zisa::int_t k = 1;
    auto stencil = zisa::biased_stencil(*grid, i_cell, k, n_points);
    REQUIRE(stencil.size() == n_points);
  }

  SECTION("biased_stencil; k = 2") {
    zisa::int_t k = 2;
    auto stencil = zisa::biased_stencil(*grid, i_cell, k, n_points);
    REQUIRE(stencil.size() == n_points);
  }
}

TEST_CASE("Stencil API", "[stencil]") {
  // FIXME remove argument
  auto grid = zisa::load_gmsh("grids/dbg.msh", 0);

  zisa::int_t i_cell = 6;

  auto params = zisa::StencilParams(3, "c", 2.0);
  auto l2g = std::vector<zisa::int_t>();

  SECTION("central_stencil") {
    // this is how one should generate a stencil
    auto approx = zisa::Stencil(l2g, grid, i_cell, params);
    REQUIRE(approx.size() > 0);

    // this is internal API, just a vector of int_t.
    auto exact = zisa::central_stencil(*grid, i_cell, approx.size());

    REQUIRE(approx.size() == exact.size());
    for (zisa::int_t i = 0; i < approx.size(); ++i) {
      REQUIRE(approx.global(i) == exact[i]);
    }
  }

  SECTION("compatibility with std::vector") {

    // FIXME remove argument
    auto grid = zisa::load_gmsh("grids/dbg.msh", 0);
    zisa::int_t i_cell = 6;
    zisa::int_t n_stencils = 3;

    auto params = zisa::StencilParams(3, "c", 2.0);

    auto l2g = std::vector<zisa::int_t>();
    auto stencils = std::vector<zisa::Stencil>();

    auto exact = zisa::Stencil(l2g, grid, i_cell, params);

    for (zisa::int_t i = 0; i < n_stencils; ++i) {
      stencils.emplace_back(l2g, grid, i_cell, params);
    }

    REQUIRE(stencils.size() == n_stencils);
    for (zisa::int_t i = 0; i < n_stencils; ++i) {
      const auto &approx = stencils[i];

      REQUIRE(approx.size() == exact.size());
      REQUIRE(approx == exact);
    }
  }
}

TEST_CASE("deduce_max_order", "[stencil]") {
  double factor = 2.0;

  int n_dims = 2;

  REQUIRE(zisa::deduce_max_order(1, factor, n_dims) == 1);
  REQUIRE(zisa::deduce_max_order(2, factor, n_dims) == 1);

  REQUIRE(zisa::deduce_max_order(4, factor, n_dims) == 2);
  REQUIRE(zisa::deduce_max_order(9, factor, n_dims) == 2);

  REQUIRE(zisa::deduce_max_order(10, factor, n_dims) == 3);
  REQUIRE(zisa::deduce_max_order(17, factor, n_dims) == 3);

  REQUIRE(zisa::deduce_max_order(18, factor, n_dims) == 4);
  REQUIRE(zisa::deduce_max_order(27, factor, n_dims) == 4);

  REQUIRE(zisa::deduce_max_order(28, factor, n_dims) == 5);
}
