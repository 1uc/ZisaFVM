#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/reconstruction/hybrid_weno.hpp>

TEST_CASE("CWENO_AO API", "[weno_ao][math]") {

  SECTION("compatibility with std::vector") {
    SECTION("push_back") {

      auto grid = zisa::load_grid("grids/small.msh", 1);
      auto params
          = zisa::HybridWENOParams({{{1}, {"c"}, {2.0}}, {1.0}, 1e-6, 4});

      auto rc = std::vector<zisa::CWENO_AO>();
      for (const auto &[i, cell] : cells(*grid)) {
        rc.emplace_back(grid, i, params);
      }

      REQUIRE(rc.size() == grid->n_cells);

      for (decltype(rc.size()) i = 0; i < rc.size(); ++i) {
        const auto &approx = rc[i];
        auto exact = zisa::CWENO_AO(grid, i, params);

        REQUIRE(approx == exact);
      }
    }
  }
}

TEST_CASE("CWENO; reconstruct smooth", "[weno_ao][math]") {
  auto grid_names
      = std::vector<std::string>{"grids/convergence/unit_square_1.msh",
                                 "grids/convergence/unit_square_2.msh"};

  using interval_t = std::tuple<double, double>;
  auto cases = std::vector<std::tuple<interval_t, zisa::HybridWENOParams>>{};

  double eps = 1e-10;
  double s = 4.0;

  cases.push_back({{2.8, 3.35},
                    {{{3, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
                     {100.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  // CWENO is expected to be high-order even for small weights of the central
  // stencil.
  cases.push_back({{3.8, 5.4},
                   {{{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
                    {10.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  cases.push_back({{3.8, 4.7},
                   {{{4, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
                    {100.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  cases.push_back({{4.8, 5.7},
                   {{{5, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
                    {100.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  for (auto &[expected_rate, params] : cases) {
    zisa::test_hybrid_weno_convergence<zisa::CWENO_AO>(
        grid_names, expected_rate, params);

    zisa::test_hybrid_weno_stability<zisa::CWENO_AO>(grid_names, params);
  }
}

TEST_CASE("CWENO; reconstruct smooth 3D", "[weno_ao][3d][math]") {
  auto grid_names = std::vector<std::string>{
      "grids/convergence/unit_cube_0.msh",
      //      "grids/convergence/unit_cube_1.msh"
  };

  using interval_t = std::tuple<double, double>;
  auto cases = std::vector<std::tuple<interval_t, zisa::HybridWENOParams>>{};

  double eps = 1e-10;
  double s = 4.0;

  cases.push_back(
      {{2.8, 3.35},
       {{{3, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"}, {2.5, 2.5, 2.5, 2.5, 2.5}},
        {100.0, 1.0, 1.0, 1.0, 1.0},
        eps,
        s}});
  cases.push_back({{2.8, 3.35}, {{{3}, {"c"}, {2.5}}, {1.0}, eps, s}});
  cases.push_back({{1.8, 2.2}, {{{2}, {"c"}, {3.0}}, {1.0}, eps, s}});

  for (auto &[expected_rate, params] : cases) {
    zisa::test_hybrid_weno_convergence<zisa::CWENO_AO>(
        grid_names, expected_rate, params);

    //    zisa::test_hybrid_weno_stability<zisa::CWENO_AO>(grid_names, params);
  }
}
