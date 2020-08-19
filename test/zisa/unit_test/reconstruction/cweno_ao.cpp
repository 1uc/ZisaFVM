#include <numeric>

#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>
#include <zisa/unit_test/reconstruction/hybrid_weno.hpp>

TEST_CASE("CWENO_AO API", "[weno_ao][math]") {
  SECTION("compatibility with std::vector") {
    SECTION("push_back") {

      auto grid = zisa::load_grid(zisa::TestGridFactory::small(), 1);
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

auto cweno_ao_2d_grid_names() {
  return std::vector<std::string>{
      zisa::TestGridFactory::unit_square_with_halo(1),
      zisa::TestGridFactory::unit_square_with_halo(2)};
}

auto cweno_ao_2d_cases() {

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
  cases.push_back({{3.8, 5.5},
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

  return cases;
}

TEST_CASE("CWENO; reconstruct smooth (stencil)",
          "[runme][weno_ao][math][2d][stencil]") {
  auto grid_names = cweno_ao_2d_grid_names();
  auto cases = cweno_ao_2d_cases();

  for (auto &[expected_rate, params] : cases) {
    zisa::test_hybrid_weno_valid_stencil<zisa::CWENO_AO>(
        grid_names, expected_rate, params, 2);
  }
}

TEST_CASE("CWENO; reconstruct smooth (rate)", "[weno_ao][math][2d][rate]") {
  auto grid_names = cweno_ao_2d_grid_names();
  auto cases = cweno_ao_2d_cases();

  for (auto &[expected_rate, params] : cases) {
    zisa::test_hybrid_weno_convergence<zisa::CWENO_AO>(
        grid_names, expected_rate, params, 2);
  }
}

TEST_CASE("CWENO; reconstruct smooth (stability)",
          "[weno_ao][math][2d][stability]") {
  auto grid_names = cweno_ao_2d_grid_names();
  auto cases = cweno_ao_2d_cases();

  for (auto &[expected_rate, params] : cases) {
    zisa::test_hybrid_weno_stability<zisa::CWENO_AO>(grid_names, params, 2);
  }
}

std::vector<
    std::tuple<bool, std::tuple<double, double>, zisa::HybridWENOParams>>
cweno_3d_cases() {
  using interval_t = std::tuple<double, double>;
  auto cases
      = std::vector<std::tuple<bool, interval_t, zisa::HybridWENOParams>>{};

  double eps = 1e-10;
  double s = 4.0;
  double o2 = 2.5;

  cases.push_back({false, {1.8, 2.2}, {{{2}, {"c"}, {4.0}}, {1.0}, eps, s}});
  cases.push_back({false, {2.7, 3.35}, {{{3}, {"c"}, {2.0}}, {1.0}, eps, s}});
  cases.push_back({false, {3.5, 4.35}, {{{4}, {"c"}, {3.0}}, {1.0}, eps, s}});

  cases.push_back(
      {true,
       {1.7, 2.2},
       {{{2, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"}, {4.0, o2, o2, o2, o2}},
        {100.0, 1.0, 1.0, 1.0, 1.0},
        eps,
        s}});

  cases.push_back(
      {true,
       {2.7, 3.1},
       {{{3, 2, 2, 2, 2}, {"c", "b", "b", "b", "b"}, {3.0, o2, o2, o2, o2}},
        {100.0, 1.0, 1.0, 1.0, 1.0},
        eps,
        s}});

  cases.push_back({true,
                   {3.7, 4.5},
                   {{{4, 2, 2, 2, 2, 2},
                     {"c", "c", "b", "b", "b", "b"},
                     {4.0, 4.0, o2, o2, o2, o2}},
                    {100.0, 10.0, 1.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  cases.push_back({true,
                   {3.7, 4.5},
                   {{{4, 3, 3, 3, 3, 3},
                     {"c", "c", "b", "b", "b", "b"},
                     {4.0, 4.0, o2, o2, o2, o2}},
                    {100.0, 10.0, 1.0, 1.0, 1.0, 1.0},
                    eps,
                    s}});

  return cases;
}

auto cweno_3d_grid_names() {
  return std::vector<std::string>{
      zisa::TestGridFactory::unit_cube_with_halo(0),
      zisa::TestGridFactory::unit_cube_with_halo(1)};
}

TEST_CASE("CWENO; reconstruct smooth 3D (matrices)",
          "[weno_ao][3d][math][matrices]") {
  auto grid_names = cweno_3d_grid_names();
  auto cases = cweno_3d_cases();

  for (auto &[_, expected_rate, params] : cases) {
    zisa::test_hybrid_weno_matrices<zisa::CWENO_AO>(grid_names, params, 3);
  }
}

TEST_CASE("CWENO; reconstruct smooth 3D (stencil)",
          "[weno_ao][3d][math][stencil]") {
  auto grid_names = cweno_3d_grid_names();
  auto cases = cweno_3d_cases();

  for (auto &[_, expected_rate, params] : cases) {
    zisa::test_hybrid_weno_valid_stencil<zisa::CWENO_AO>(
        grid_names, expected_rate, params, 3);
  }
}

TEST_CASE("CWENO; reconstruct smooth 3D (rate)", "[weno_ao][3d][math][rate]") {
  auto grid_names = cweno_3d_grid_names();
  auto cases = cweno_3d_cases();

  for (auto &[_, expected_rate, params] : cases) {
    zisa::test_hybrid_weno_convergence<zisa::CWENO_AO>(
        grid_names, expected_rate, params, 3);
  }
}

TEST_CASE("CWENO; reconstruct smooth 3D (stability)",
          "[weno_ao][3d][math][stability]") {
  auto grid_names = cweno_3d_grid_names();
  auto cases = cweno_3d_cases();
  for (auto &[is_stable, expected_rate, params] : cases) {
    if (is_stable) {
      zisa::test_hybrid_weno_stability<zisa::CWENO_AO>(grid_names, params, 3);
    }
  }
}
