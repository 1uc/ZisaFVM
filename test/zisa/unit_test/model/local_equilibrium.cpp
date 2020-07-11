#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/cell_factory.hpp>
#include <zisa/math/triangle.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/model/isentropic_equilibrium.hpp>
#include <zisa/model/janka_eos.hpp>
#include <zisa/model/local_equilibrium.hpp>
#include <zisa/model/polytrope.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("LocalEquilibrium", "[equilibrium]") {
  using euler_t = zisa::Euler;
  auto eos = std::make_shared<zisa::IdealGasEOS>(1.2, 0.9);
  auto gravity = std::make_shared<zisa::ConstantGravityRadial>(0.9);
  auto euler = std::make_shared<euler_t>();

  zisa::int_t quad_deg = 2;
  auto eq = zisa::IsentropicEquilibrium(eos, gravity);

  auto tri_ref
      = zisa::Triangle{{1.0, 1.0, 0.0}, {1.01, 1.0, 0.0}, {1.0, 1.01, 0.0}};

  auto cell_ref = make_cell(tri_ref, quad_deg);

  auto x_ref = zisa::XYZ{0.5, 0.6, 0.0};
  auto theta_ref = zisa::EnthalpyEntropy{10.0, 3.0};

  auto rhoE_eq = [&eq, &theta_ref, &x_ref](const zisa::XYZ &xy) {
    return eq.extrapolate(theta_ref, x_ref, xy);
  };

  auto rhoE_bar = average(cell_ref, rhoE_eq);

  auto eq_loc = zisa::LocalEquilibrium(eq);
  eq_loc.solve(rhoE_bar, cell_ref);

  SECTION("extrapolate to point") {
    auto xy = zisa::XYZ{1.1, 2.1, 0.0};

    auto approx = eq_loc.extrapolate(xy);
    auto exact = rhoE_eq(xy);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-9));
  }

  SECTION("extrapolate to triangle") {
    auto tri
        = zisa::Triangle{{1.2, 1.1, 0.0}, {1.21, 1.1, 0.0}, {1.2, 1.11, 0.0}};
    auto cell = make_cell(tri, quad_deg);

    auto approx = eq_loc.extrapolate(cell);
    auto exact = average(cell, rhoE_eq);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-9));
  }
}
