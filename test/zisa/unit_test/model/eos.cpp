#include <zisa/model/ideal_gas_eos.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("LocalEquilibrium", "[equilibrium]") {
  auto eos = zisa::IdealGasEOS(1.2, 0.9);

  SECTION("rhoE -> hK -> rhoE") {
    auto rhoE = zisa::RhoE{1.0, 2.0};
    auto hK = eos.enthalpy_entropy(rhoE);

    REQUIRE(zisa::almost_equal(eos.rhoE(hK), rhoE, 1e-10));
  }
}
