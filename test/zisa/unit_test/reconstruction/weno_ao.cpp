#include <numeric>

#include <zisa/testing/testing_framework.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/weno_ao.hpp>
#include <zisa/unit_test/reconstruction/hybrid_weno.hpp>

TEST_CASE("WENO_AO API", "[weno_ao][math]") {

  SECTION("compatibility with std::vector") {
    SECTION("push_back") {

      auto grid = zisa::load_gmsh("grids/small.msh");
      auto params
          = zisa::HybridWENOParams({{{1}, {"c"}, {2.0}}, {1.0}, 1e-10, 4.0});

      auto rc = std::vector<zisa::WENO_AO>();
      for (const auto &[i, tri] : triangles(*grid)) {
        rc.push_back(zisa::WENO_AO(grid, i, params));
      }

      REQUIRE(rc.size() == grid->n_cells);

      for (decltype(rc.size()) i = 0; i < rc.size(); ++i) {
        const auto &approx = rc[i];
        auto exact = zisa::WENO_AO(grid, i, params);

        REQUIRE(approx == exact);
      }
    }
  }
}

TEST_CASE("WENO_AO; reconstruct smooth", "[weno_ao][math]") {
  auto grid_names
      = std::vector<std::string>{"grids/convergence/unit_square_1.msh",
                                 "grids/convergence/unit_square_2.msh"};

  double eps = 1e-10;
  double s = 4.0;

  using interval_t = std::tuple<double, double>;
  auto cases
      = std::vector<std::tuple<interval_t, bool, zisa::HybridWENOParams>>{
          {{0.8, 1.15}, true, {{{1}, {"c"}, {2.0}}, {1.0}, eps, s}},
          {{0.8, 1.15}, true, {{{1}, {"b"}, {2.0}}, {1.0}, eps, s}},
          {{1.8, 2.2}, false, {{{2}, {"c"}, {2.0}}, {1.0}, eps, s}},
          {{1.8, 2.2}, false, {{{2}, {"b"}, {2.0}}, {1.0}, eps, s}},
          {{2.8, 3.25}, false, {{{3}, {"c"}, {2.0}}, {1.0}, eps, s}},
          {{2.8, 3.25}, false, {{{3}, {"b"}, {2.0}}, {1.0}, eps, s}},
          {{3.8, 4.4}, false, {{{4}, {"c"}, {2.0}}, {1.0}, eps, s}}};

  cases.push_back({{2.9, 3.3},
                   true,
                   {{{2, 2, 2, 3}, {"b", "b", "b", "c"}, {1.5, 1.5, 1.5, 2.0}},
                    {1.0, 1.0, 1.0, 100.0},
                    eps,
                    s}});

  // The reason for the second order convergence is the linear weights. They
  // allow too much pollution from the second order stencils.
  cases.push_back({{2.0, 3.8},
                   false /* empirically is not non-oscillatory */,
                   {{{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
                    {10.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  // Which can be remedied by increasing the central weight.
  cases.push_back({{3.9, 4.4},
                   false /* not non-oscillatory */,
                   {{{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
                    {1000.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  for (auto &[expected_rate, is_stable, params] : cases) {
    zisa::test_hybrid_weno_convergence<zisa::WENO_AO>(
        grid_names, expected_rate, params);

    if (is_stable) {
      zisa::test_hybrid_weno_stability<zisa::WENO_AO>(grid_names, params);
    }
  }
}
