#include <catch/catch.hpp>

#include <zisa/math/poly2d.hpp>

TEST_CASE("poly_dof", "[math][poly2d]") {
  REQUIRE(zisa::poly_dof(0) == 1);
  REQUIRE(zisa::poly_dof(1) == 3);
  REQUIRE(zisa::poly_dof(2) == 6);
  REQUIRE(zisa::poly_dof(3) == 10);
  REQUIRE(zisa::poly_dof(4) == 15);
}

TEST_CASE("poly_degree", "[math][poly2d]") {

  SECTION("right-inverse") {
    for (zisa::int_t n_coeffs = 1; n_coeffs < 25; ++n_coeffs) {
      auto deg = zisa::poly_degree(n_coeffs);

      INFO(string_format("n_coeffs = %d \n", n_coeffs));
      REQUIRE(zisa::poly_dof(deg) <= n_coeffs);
      REQUIRE(zisa::poly_dof(deg + 1) > n_coeffs);
    }
  }

  SECTION("left-inverse") {
    for (zisa::int_t deg = 0; deg < 5; ++deg) {
      INFO(string_format("deg = %d \n", deg));
      REQUIRE(zisa::poly_degree(zisa::poly_dof(deg)) == deg);
    }
  }
}

TEST_CASE("Poly2D; examples", "[math][poly2d]") {
  auto p = zisa::Poly2D<5>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                           {0.0, 0.0, 0.0, 1.0, 2.0, 3.0});

  SECTION("degree") { REQUIRE(p.degree() == 2); }
  SECTION("max_degree") { REQUIRE(p.max_degree() == 5); }

  SECTION("x == 0") {
    auto approx = p(zisa::XY{0.0, 0.0});
    auto exact = 1.0 - 4.0 - 10.0 - 18.0;

    INFO(string_format("%e != %e [%e]\n", approx, exact, approx - exact));
    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }

  SECTION("x != 0") {
    double x = -3.0, y = 2.0;

    auto approx = p(zisa::XY{x, y});
    auto exact = 1.0 + 2.0 * x + 3.0 * y + 4.0 * (x * x - 1.0)
                 + 5.0 * (x * y - 2.0) + 6.0 * (y * y - 3.0);

    INFO(string_format("%e != %e [%e]\n", approx, exact, approx - exact));
    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }

  SECTION("saxpy-like") {
    auto p = zisa::Poly2D<5>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                             {0.0, 0.0, 0.0, 1.0, 2.0, 3.0});

    auto q = zisa::Poly2D<5>({1.0, 2.0, 3.0},
                             {0.0, 0.0, 0.0});

    auto x = zisa::XY{-3.4, 2.138};

    auto pq = zisa::Poly2D<5>(0.2*p + q - 0.4*p);

    auto exact = 0.2*p(x) + q(x) - 0.4*p(x);
    auto approx = pq(x);

    INFO(string_format("%e != %e [%e]\n", approx, exact, approx - exact));
    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }
}
