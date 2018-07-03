#include <catch/catch.hpp>

#include <zisa/flux/hllc.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/ideal_gas_eos.hpp>

TEST_CASE("HLLC; consistency") {

  using eos_t = zisa::IdealGasEOS;
  using gravity_t = zisa::ConstantGravityRadial;
  using euler_t = zisa::Euler<eos_t, gravity_t>;

  auto eos = eos_t{1.6, 1.0};
  auto gravity = gravity_t{0.9};

  auto euler = euler_t(0.0, eos, gravity);
  auto u = zisa::euler_var_t{1.0, -0.2, 0.3, 0.8, 12.0};
  auto p = eos.pressure(u);

  auto nf = zisa::HLLCBatten<euler_t>::flux(euler, u, u);
  auto pf = euler.flux(u, p);

  REQUIRE(zisa::almost_equal(nf, pf, 1e-12));
}
