#include <zisa/math/triangle.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/local_equilibrium.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("LocalEquilibrium", "[equilibrium]") {
  auto eos = zisa::IdealGasEOS(1.2, 0.9);
  auto gravity = zisa::ConstantGravityRadial(0.9);

  zisa::int_t quad_deg = 2;
  auto eq = zisa::IsentropicEquilibrium(eos, gravity, quad_deg);

  auto tri_ref = zisa::Triangle{{1.0, 1.0}, {1.01, 1.0}, {1.0, 1.01}};
  auto x_ref = zisa::XY{0.5, 0.6};
  auto theta_ref = zisa::EnthalpyEntropy{10.0, 3.0};

  auto rhoE_eq = [&eq, &theta_ref, &x_ref](const zisa::XY &xy) {
    return extrapolate(eq, theta_ref, x_ref, xy);
  };

  auto vol = volume(tri_ref);
  auto rhoE_bar = average(rhoE_eq, tri_ref, quad_deg);

  auto eq_loc = zisa::LocalEquilibrium(eq, tri_ref);
  eq_loc.solve(rhoE_bar);

  SECTION("extrapolate to point") {
    auto xy = zisa::XY{1.1, 2.1};

    auto approx = eq_loc.extrapolate(xy);
    auto exact = rhoE_eq(xy);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-9));
  }

  SECTION("extrapolate to triangle") {
    auto tri = zisa::Triangle{{1.2, 1.1}, {1.21, 1.1}, {1.2, 1.11}};
    auto vol = zisa::volume(tri);

    auto approx = eq_loc.extrapolate(tri);
    auto exact = average(rhoE_eq, tri, quad_deg);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-9));
  }
}
