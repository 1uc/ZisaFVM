// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/reconstruction/stencil_family.hpp>

#include <zisa/math/poly2d.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

void check_uniqueness(const zisa::StencilFamily &stencils, zisa::int_t i_cell) {
  const auto &l2g = stencils.local2global();

  REQUIRE(l2g[0] == i_cell);
  for (zisa::int_t i = 1; i < l2g.size(); ++i) {
    REQUIRE(l2g[i] != i_cell);
  }
}

TEST_CASE("StencilFamily", "[weno_ao]") {

  SECTION("initialization") {
    auto grid = zisa::load_grid(zisa::TestGridFactory::unit_cube(0), 3);

    SECTION("is_unique") {
      for (zisa::int_t i_cell = 0; i_cell < grid->n_cells; ++i_cell) {
        auto biases = std::vector<std::string>{"c", "b"};
        for (auto &&b : biases) {
          auto o1_stencils
              = zisa::StencilFamily(*grid, i_cell, {{1}, {b}, {2.0}});
          INFO(string_format("[% 4d] bias = %s", i_cell, b.c_str()));
          check_uniqueness(o1_stencils, i_cell);

          auto o2_stencils
              = zisa::StencilFamily(*grid, i_cell, {{2}, {b}, {2.0}});
          INFO(string_format("[% 4d] bias = %s", i_cell, b.c_str()));
          check_uniqueness(o2_stencils, i_cell);

          auto o3_stencils
              = zisa::StencilFamily(*grid, i_cell, {{3}, {b}, {2.0}});
          INFO(string_format("[% 4d] bias = %s", i_cell, b.c_str()));
          check_uniqueness(o3_stencils, i_cell);
        }
      }
    }

    zisa::int_t i_cell = 20;
    SECTION("single stencil") {
      auto biases = std::vector<std::string>{"c", "b"};
      for (auto &&b : biases) {
        auto o1_stencils
            = zisa::StencilFamily(*grid, i_cell, {{1}, {b}, {2.0}});
        INFO(string_format("bias = %s", b.c_str()));
        REQUIRE(o1_stencils.order() == 1);

        auto o2_stencils
            = zisa::StencilFamily(*grid, i_cell, {{2}, {b}, {2.0}});
        INFO(string_format("bias = %s", b.c_str()));
        REQUIRE(o2_stencils.order() == 2);

        auto o3_stencils
            = zisa::StencilFamily(*grid, i_cell, {{3}, {b}, {2.0}});
        INFO(string_format("bias = %s", b.c_str()));
        REQUIRE(o3_stencils.order() == 3);
      }
    }

    SECTION("mixed stencil") {
      auto o311_stencils = zisa::StencilFamily(
          *grid, i_cell, {{3, 1, 1}, {"c", "b", "b"}, {2.0, 1.5, 1.5}});

      REQUIRE(o311_stencils.order() == 3);
    }
  }

  SECTION("compatibility with std::vector") {
    SECTION("push_back") {
      auto grid = zisa::load_grid(zisa::TestGridFactory::small());
      auto params = zisa::StencilFamilyParams({{1}, {"c"}, {2.0}});

      auto stencils = std::vector<zisa::StencilFamily>();
      for (const auto &[i, cell] : cells(*grid)) {
        stencils.emplace_back(*grid, i, params);
      }

      REQUIRE(stencils.size() == grid->n_cells);
    }
  }
}

TEST_CASE("StencilFamily, 3D", "[stencil][3d]") {
  auto grid = zisa::load_grid(zisa::TestGridFactory::unit_cube(1), 3);
  zisa::int_t i_cell = 90;
  auto n_cells = grid->n_cells;

  auto cell_indices = std::vector<zisa::int_t>{0, 20, 90};

  SECTION("mixed stencil") {
    for (auto i_cell : cell_indices) {
      auto stencils = zisa::StencilFamily(*grid,
                                          i_cell,
                                          {{3, 2, 2, 2, 2},
                                           {"c", "b", "b", "b", "b"},
                                           {2.0, 1.5, 1.5, 1.5, 1.5}});

      for (zisa::int_t k = 0; k < 5; ++k) {
        const auto &s = stencils[k];

        auto n = zisa::poly_dof<3>(s.order() - 1);
        auto observed = double(s.size());
        auto expected = double(n * s.overfit_factor());

        // Poly has zero mean.
        INFO(string_format("i = %d", i_cell));
        REQUIRE(zisa::abs(observed - expected)
                <= zisa::ceil(s.overfit_factor()));
      }

      INFO(string_format("i = %d", i_cell));
      REQUIRE(stencils.order() == 3);
    }
  }
}

TEST_CASE("StencilFamily; real-world issue", "[2d][stencil]") {
  auto grid = zisa::load_grid(zisa::TestGridFactory::unit_cube_with_halo(1), 3);

  zisa::int_t i = 24;

  auto sf = zisa::StencilFamily(
      *grid, i, {{3, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}});

  REQUIRE(sf[0].order() == 3);
  REQUIRE(sf[1].order() == 2);
  REQUIRE(sf[2].order() == 2);
  REQUIRE(sf[3].order() == 2);
}
