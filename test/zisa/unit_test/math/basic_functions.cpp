#include "basic_functions.hpp"
#include <zisa/math/basic_functions.hpp>

#include <cassert>
#include <catch/catch.hpp>
#include <cmath>

std::vector<double> convergence_rates(const std::vector<double> &dx,
                                      const std::vector<double> &e) {

  assert(dx.size() == e.size());
  auto r = std::vector<double>(e.size() - 1);

  for (std::size_t i = 0; i < r.size(); ++i) {
    r[i] = (std::log(e[i + 1]) - std::log(e[i]))
           / (std::log(dx[i + 1]) - std::log(dx[i]));
  }

  return r;
}

TEST_CASE("factorial") {
  REQUIRE(zisa::factorial(0) == 1);
  REQUIRE(zisa::factorial(1) == 1);
  REQUIRE(zisa::factorial(2) == 2);
  REQUIRE(zisa::factorial(3) == 6);
  REQUIRE(zisa::factorial(4) == 24);
}

TEST_CASE("Gamma") {
  for (int k = 1; k < 6; ++k) {
    REQUIRE(zisa::Gamma(k) == zisa::factorial(k - 1));
  }
}
