#include <zisa/math/cone.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("Cone", "[cone][math]") {

  SECTION("cone A") {
    auto c = zisa::Cone({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0});

    REQUIRE(c.is_inside({0.5, 0.5, 0.0}));
    REQUIRE(!c.is_inside({-0.5, 0.5, 0.0}));
    REQUIRE(!c.is_inside({0.5, -0.5, 0.0}));
    REQUIRE(!c.is_inside({-0.5, -0.5, 0.0}));
  }

  SECTION("cone B") {
    auto c = zisa::Cone({-0.0, -2.0, 0.0}, {1.0, -2.0, 0.0}, {1.0, 1.0, 0.0});

    REQUIRE(c.is_inside({0.5, 0.5, 0.0}));
    REQUIRE(c.is_inside({0.5, 0.5, 0.1}));
    REQUIRE(c.is_inside({-0.5, 0.5, 0.0}));
    REQUIRE(c.is_inside({-0.5, 0.5, 0.1}));

    REQUIRE(!c.is_inside({3.5, -0.5, 0.0}));
    REQUIRE(!c.is_inside({3.5, -0.5, 1.0}));
    REQUIRE(!c.is_inside({0.5, -4.5, 0.0}));
    REQUIRE(!c.is_inside({0.5, -4.5, 10.0}));
  }
}
