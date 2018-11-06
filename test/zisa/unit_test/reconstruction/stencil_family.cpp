#include <catch/catch.hpp>

#include <zisa/reconstruction/stencil_family.hpp>

TEST_CASE("StencilFamily", "[weno_ao]") {

  SECTION("initialization") {
    auto grid = zisa::load_gmsh("grids/small.msh");
    zisa::int_t i_cell = 20;

    SECTION("single stencil") {
      auto biases = std::vector<std::string>{"c", "b"};
      for (auto &&b : biases) {
        auto o1_stencils = zisa::StencilFamily(grid, i_cell, {{1}, {b}, {2.0}});
        INFO(string_format("bias = %s", b.c_str()));
        REQUIRE(o1_stencils.order() == 1);

        auto o2_stencils = zisa::StencilFamily(grid, i_cell, {{2}, {b}, {2.0}});
        INFO(string_format("bias = %s", b.c_str()));
        REQUIRE(o2_stencils.order() == 2);

        auto o3_stencils = zisa::StencilFamily(grid, i_cell, {{3}, {b}, {2.0}});
        INFO(string_format("bias = %s", b.c_str()));
        REQUIRE(o3_stencils.order() == 3);
      }
    }

    SECTION("mixed stencil") {
      auto o311_stencils = zisa::StencilFamily(
          grid, i_cell, {{3, 1, 1}, {"c", "b", "b"}, {2.0, 1.5, 1.5}});

      REQUIRE(o311_stencils.order() == 3);
    }
  }

  SECTION("compatibility with std::vector") {
    SECTION("push_back") {

      auto grid = zisa::load_gmsh("grids/small.msh");
      auto params = zisa::StencilFamilyParams({{1}, {"c"}, {2.0}});

      auto stencils = std::vector<zisa::StencilFamily>();
      for (const auto &[i, tri] : triangles(*grid)) {
        stencils.push_back(zisa::StencilFamily(grid, i, params));
      }

      REQUIRE(stencils.size() == grid->n_cells);

      // REQUIRE(stencils.size() == n_stencils);
      // for (zisa::int_t i = 0; i < n_stencils; ++i) {
      //   const auto &approx = stencils[i];

      //   REQUIRE(approx.size() == exact.size());
      //   REQUIRE(approx == exact);
      // }
    }


  }
}