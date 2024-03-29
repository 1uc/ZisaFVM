// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <random>

#include <zisa/testing/testing_framework.hpp>

#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

TEST_CASE("Wellbalanced RC; small perturbations", "[wb][math]") {

  auto quad_deg = zisa::int_t(2);
  auto grid = zisa::load_grid(zisa::TestGridFactory::polytrope(), quad_deg);
  auto n_cells = grid->n_cells;
  auto n_cvars = zisa::int_t(5);
  auto n_avars = zisa::int_t(0);

  auto dims = zisa::AllVariablesDimensions{n_cells, n_cvars, n_avars};
  auto all_variables = std::make_shared<zisa::AllVariables>(dims);

  using eos_t = zisa::IdealGasEOS;
  using gravity_t = zisa::PolytropeGravityRadial;
  using euler_t = zisa::Euler;

  auto local_eos = zisa::LocalEOSState<eos_t>(2.0, 1.0);
  auto eos = local_eos(0);
  auto gravity = std::make_shared<gravity_t>();
  auto euler = std::make_shared<euler_t>();

  double width = 0.05;
  double amp = 0.0;
  double rand_amp = 1e-8;
  double atol = 3.3 * rand_amp;

  auto rd = std::random_device();
  auto gen = std::mt19937(rd());
  auto dist = std::uniform_int_distribution(-1, 1);

  auto ic_ = zisa::PolytropeIC(eos, gravity);
  auto ic = [&eos, &ic_, amp, width](const auto &x) {
    // Avoid implicit conversion.
    zisa::RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;
    double p = p_eq * (1 + amp * zisa::exp(-zisa::pow<2>(r / width)));

    return eos->cvars(zisa::RhoP{rho_eq, p});
  };

  auto &u0 = all_variables->cvars;
  for (auto &&[i, cell] : zisa::cells(*grid)) {
    u0(i) = zisa::average(cell, ic);
    u0(i, 4) += dist(gen) * rand_amp;
  }

  double eps = 1e-10;
  double s = 4.0;

  auto weno_params = zisa::HybridWENOParams(
      {{2, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
      {100.0, 1.0, 1.0, 1.0},
      eps,
      s);

  using eq_t = zisa::IsentropicEquilibrium<eos_t, gravity_t>;
  auto eq = eq_t(eos, gravity);

  using rc_t = zisa::CWENO_AO;
  using scaling_t = zisa::UnityScaling;
  auto scaling = scaling_t{};
  auto local_rc_params = zisa::LocalRCParams{1, -1.0};
  auto rc = zisa::
      make_reconstruction_array<eq_t, rc_t, scaling_t, eos_t, gravity_t>(
          grid, weno_params, local_eos, gravity, local_rc_params);

  auto grc = zisa::EulerGlobalReconstruction(weno_params, rc);

  grc.compute(*all_variables);

  for (auto &&[i, cell] : zisa::cells(*grid)) {
    for (const auto &x : cell.qr.points) {
      double dE = zisa::norm(grc(i, x) - ic(x));

      INFO(string_format("[%d] %.3e >= %.3e", i, dE, atol));
      REQUIRE(dE < atol);
    }
  }
}
