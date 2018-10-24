#include <catch/catch.hpp>

#include <zisa/reconstruction/stencil_family.hpp>

TEST_CASE("StencilFamily", "[weno_ao]") {

  SECTION("initialization") {
    auto grid = zisa::load_gmsh("grids/small.msh");
    zisa::int_t i_cell = 20;

    SECTION("single stencil") {
      auto o1_stencils = zisa::StencilFamily(grid, i_cell, {{1}, {"c"}, {2.0}});
      REQUIRE(o1_stencils.order() == 1);

      auto o2_stencils = zisa::StencilFamily(grid, i_cell, {{2}, {"c"}, {2.0}});
      REQUIRE(o2_stencils.order() == 2);

      auto o3_stencils = zisa::StencilFamily(grid, i_cell, {{3}, {"c"}, {2.0}});
      REQUIRE(o3_stencils.order() == 3);
    }

    SECTION("mixed stencil") {
      auto o311_stencils = zisa::StencilFamily(
          grid, i_cell, {{3, 1, 1}, {"c", "b", "b"}, {2.0, 1.5, 1.5}});

      REQUIRE(o311_stencils.order() == 3);
    }
  }
}
