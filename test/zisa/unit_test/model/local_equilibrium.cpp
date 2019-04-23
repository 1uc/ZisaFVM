#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/janka_eos.hpp>
#include <zisa/model/local_equilibrium.hpp>
#include <zisa/model/polytrope.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("LocalEquilibrium", "[equilibrium]") {
  using euler_t = zisa::Euler<zisa::IdealGasEOS, zisa::ConstantGravityRadial>;
  auto eos = zisa::IdealGasEOS(1.2, 0.9);
  auto gravity = zisa::ConstantGravityRadial(0.9);
  auto euler = std::make_shared<euler_t>(eos, gravity);

  zisa::int_t quad_deg = 2;
  auto eq = zisa::IsentropicEquilibrium(euler, quad_deg);

  auto tri_ref
      = zisa::Triangle{{1.0, 1.0, 0.0}, {1.01, 1.0, 0.0}, {1.0, 1.01, 0.0}};
  auto x_ref = zisa::XYZ{0.5, 0.6, 0.0};
  auto theta_ref = zisa::EnthalpyEntropy{10.0, 3.0};

  auto rhoE_eq = [&eq, &theta_ref, &x_ref](const zisa::XYZ &xy) {
    return eq.extrapolate(theta_ref, x_ref, xy);
  };

  auto rhoE_bar = average(rhoE_eq, tri_ref, quad_deg);

  auto eq_loc = zisa::LocalEquilibrium(eq);
  eq_loc.solve(rhoE_bar, tri_ref);

  SECTION("extrapolate to point") {
    auto xy = zisa::XYZ{1.1, 2.1, 0.0};

    auto approx = eq_loc.extrapolate(xy);
    auto exact = rhoE_eq(xy);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-9));
  }

  SECTION("extrapolate to triangle") {
    auto tri
        = zisa::Triangle{{1.2, 1.1, 0.0}, {1.21, 1.1, 0.0}, {1.2, 1.11, 0.0}};

    auto approx = eq_loc.extrapolate(tri);
    auto exact = average(rhoE_eq, tri, quad_deg);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-9));
  }
}
