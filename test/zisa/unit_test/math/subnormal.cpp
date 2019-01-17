#include <zisa/testing/testing_framework.hpp>

#include <type_traits>

TEST_CASE("Inverse of normal finite numbers.") {

  SECTION("float") {
    float x_min = std::numeric_limits<float>::min();
    float x_max = std::numeric_limits<float>::max();

    REQUIRE(std::isnormal(1.0f / x_min));

    INFO(1.0f / x_max);
    REQUIRE(std::isfinite(1.0f / x_max));  // but subnormal
  }

  SECTION("double") {
    double x_min = std::numeric_limits<double>::min();
    double x_max = std::numeric_limits<double>::max();

    REQUIRE(std::isnormal(1.0 / x_min));

    INFO(1.0 / x_max);
    REQUIRE(std::isfinite(1.0 / x_max));  // but subnormal
  }
}
