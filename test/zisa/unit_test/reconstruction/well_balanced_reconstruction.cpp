#include <random>

#include <zisa/testing/testing_framework.hpp>

#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

TEST_CASE("Wellbalanced RC; small perturbations", "[wb][math]") {

  auto grid = zisa::load_gmsh("grids/unit_tests/polytrope.msh");
  auto n_cells = grid->n_cells;
  auto n_cvars = zisa::int_t(5);
  auto n_avars = zisa::int_t(0);

  auto dims = zisa::AllVariablesDimensions{n_cells, n_cvars, n_avars};
  auto all_variables = std::make_shared<zisa::AllVariables>(dims);

  auto quad_deg = zisa::int_t(2);
  auto qr = zisa::cached_triangular_quadrature_rule(quad_deg);

  using eos_t = zisa::IdealGasEOS;
  using gravity_t = zisa::PolytropeGravityRadial;

  auto euler = zisa::Euler<eos_t, gravity_t>({2.0, 1.0}, {});
  const auto &eos = euler.eos;

  double width = 0.05;
  double amp = 0.0;
  double rand_amp = 1e-8;
  double atol = 3.3 * rand_amp;

  auto rd = std::random_device();
  auto gen = std::mt19937(rd());
  auto dist = std::uniform_int_distribution(-1, 1);

  auto ic_ = zisa::PolytropeIC(euler);
  auto ic = [&eos, &ic_, amp, width](const auto &x) {
    // Avoid implicit conversion.
    zisa::RhoP rhoP = ic_(x);

    double r = zisa::norm(x);
    auto &[rho_eq, p_eq] = rhoP;
    double p = p_eq * (1 + amp * zisa::exp(-zisa::pow<2>(r / width)));

    return eos.cvars(zisa::RhoP{rho_eq, p});
  };

  auto &u0 = all_variables->cvars;
  for (auto &&[i, tri] : zisa::triangles(*grid)) {
    u0(i) = zisa::average(qr, ic, tri);
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
  auto eq = eq_t{euler.eos, euler.gravity, quad_deg};

  auto grc = zisa::EulerGlobalReconstruction<eq_t, zisa::CWENO_AO>(
      grid, weno_params, eq);

  grc.compute(*all_variables);

  for (auto &&[i, tri] : zisa::triangles(*grid)) {
    for (const auto &chi : qr.points) {
      auto x = zisa::coord(tri, chi);
      double dE = zisa::norm(grc(i, x) - ic(x));

      INFO(string_format("[%d] %.3e >= %.3e", i, dE, atol));
      REQUIRE(dE < atol);
    }
  }
}
