#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/poly2d.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/math/triangle.hpp>

TEST_CASE("poly_dof", "[math][poly2d]") {
  REQUIRE(zisa::poly_dof(0) == 1);
  REQUIRE(zisa::poly_dof(1) == 3);
  REQUIRE(zisa::poly_dof(2) == 6);
  REQUIRE(zisa::poly_dof(3) == 10);
  REQUIRE(zisa::poly_dof(4) == 15);
}

TEST_CASE("poly_index", "[math][poly2d]") {

  int i = 0;

  REQUIRE(zisa::poly_index(0, 0) == i++);

  REQUIRE(zisa::poly_index(1, 0) == i++);
  REQUIRE(zisa::poly_index(0, 1) == i++);

  REQUIRE(zisa::poly_index(2, 0) == i++);
  REQUIRE(zisa::poly_index(1, 1) == i++);
  REQUIRE(zisa::poly_index(0, 2) == i++);

  REQUIRE(zisa::poly_index(3, 0) == i++);
  REQUIRE(zisa::poly_index(2, 1) == i++);
  REQUIRE(zisa::poly_index(1, 2) == i++);
  REQUIRE(zisa::poly_index(0, 3) == i++);
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
  auto p = zisa::Poly2D<5, 1>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                              {0.0, 0.0, 0.0, 1.0, 2.0, 3.0},
                              zisa::XYZ::zeros(), 1.0);

  auto q = zisa::Poly2D<5, 1>({1.0, 2.0, 3.0}, {0.0, 0.0, 0.0},
                              zisa::XYZ::zeros(), 1.0);

  SECTION("degree") { REQUIRE(p.degree() == 2); }
  SECTION("max_degree") { REQUIRE(p.max_degree() == 5); }
  SECTION("n_vars") { REQUIRE(p.n_vars() == 1); }

  SECTION("x == 0") {
    auto approx = p(zisa::XYZ::zeros());
    auto exact = 1.0 - 4.0 - 10.0 - 18.0;

    INFO(string_format("%e != %e [%e]\n", approx[0], exact, approx - exact));
    REQUIRE(zisa::almost_equal(approx[0], exact, 1e-14));
  }

  SECTION("x != 0") {
    double x = -3.0, y = 2.0, z = 0.0;

    auto approx = p(zisa::XYZ{x, y, z});
    auto exact = 1.0 + 2.0 * x + 3.0 * y + 4.0 * (x * x - 1.0) +
                 5.0 * (x * y - 2.0) + 6.0 * (y * y - 3.0);

    INFO(string_format("%e != %e [%e]\n", approx[0], exact, approx - exact));
    REQUIRE(zisa::almost_equal(approx[0], exact, 1e-14));
  }

  SECTION("saxpy-like") {
    auto x = zisa::XYZ{-3.4, 2.138, 0.0};

    auto pq = zisa::Poly2D<5, 1>(0.2 * p + q - 0.4 * p);

    auto exact = zisa::Cartesian<1>{0.2 * p(x) + q(x) - 0.4 * p(x)};
    auto approx = pq(x);

    INFO(string_format("%e != %e [%e]\n", approx[0], exact[0], approx - exact));
    REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
  }

  SECTION("assignment") {
    auto x = zisa::XYZ{-3.4, 2.138, 0.0};

    auto tmp = zisa::Poly2D<5, 1>(p);

    SECTION("+=") {
      tmp += q;

      auto exact = zisa::Cartesian<1>{p(x) + q(x)};
      auto approx = tmp(x);

      INFO(string_format("%e != %e [%e]\n", approx[0], exact[0],
                         approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }

    SECTION("+=") {
      tmp -= q;

      auto exact = zisa::Cartesian<1>{p(x) - q(x)};
      auto approx = zisa::Cartesian<1>{tmp(x)};

      INFO(string_format("%e != %e [%e]\n", approx[0], exact[0],
                         approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }

    SECTION("*=") {
      tmp *= 2.0;

      auto exact = zisa::Cartesian<1>{2.0 * p(x)};
      auto approx = zisa::Cartesian<1>{tmp(x)};

      INFO(string_format("%e != %e [%e]\n", approx[0], exact[0],
                         approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }

    SECTION("/=") {
      tmp /= 2.0;

      auto exact = zisa::Cartesian<1>{0.5 * p(x)};
      auto approx = zisa::Cartesian<1>{tmp(x)};

      INFO(string_format("%e != %e [%e]\n", approx[0], exact[0],
                         approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }
  }

  SECTION("quadrature") {
    auto tri = zisa::reference_triangle();

    auto pbar = zisa::quadrature(p, tri, 4);
    REQUIRE(pbar[0] < 0.0);
  }
}

TEST_CASE("Poly2D<., 2>; examples", "[math][poly2d]") {
  constexpr int n_vars = 2;
  auto p = zisa::Poly2D<5, n_vars>(2, {0.0, 0.0, 0.0, 1.0, 2.0, 3.0},
                                   zisa::XYZ::zeros(), 1.0);

  auto p00 = zisa::Cartesian<n_vars>{1.0, -1.0};
  auto p01 = zisa::Cartesian<n_vars>{-1.0, -1.0};
  auto p10 = zisa::Cartesian<n_vars>{2.0, -1.0};
  auto p02 = zisa::Cartesian<n_vars>{1.0, -3.0};
  auto p11 = zisa::Cartesian<n_vars>{4.0, -1.0};
  auto p20 = zisa::Cartesian<n_vars>{0.2, -0.1};

  auto p_coeffs =
      std::vector<zisa::Cartesian<n_vars>>{p00, p10, p01, p20, p11, p02};

  std::copy((double *)&p_coeffs[0],
            (double *)&p_coeffs[0] + p.n_vars() * p.dof(p.degree()),
            p.coeffs_ptr());

  auto q = zisa::Poly2D<5, n_vars>(1, {0.0, 0.0, 0.0}, zisa::XYZ::zeros(), 1.0);

  auto q00 = zisa::Cartesian<n_vars>{1.0, -1.0};
  auto q01 = zisa::Cartesian<n_vars>{-1.0, -1.0};
  auto q10 = zisa::Cartesian<n_vars>{2.0, -1.0};

  auto q_coeffs = std::vector<zisa::Cartesian<n_vars>>{q00, q10, q01};
  std::copy((double *)&q_coeffs[0],
            (double *)&q_coeffs[0] + q.n_vars() * q.dof(q.degree()),
            q.coeffs_ptr());

  SECTION("degree") { REQUIRE(p.degree() == 2); }
  SECTION("max_degree") { REQUIRE(p.max_degree() == 5); }
  SECTION("n_vars") { REQUIRE(p.n_vars() == n_vars); }

  SECTION("p{} == 0") {
    auto p = std::make_shared<zisa::Poly2D<5, 2>>();
    REQUIRE((*p)(zisa::XYZ{2.0, -1.0, 0.0}) == zisa::Cartesian<2>{0.0, 0.0});
  }

  SECTION("p(x)") {
    double x = -3.0, y = 2.0, z = 0.0;

    auto approx = p(zisa::XYZ{x, y, z});
    auto exact =
        zisa::Cartesian<n_vars>(p00 + p10 * x + p01 * y + p20 * (x * x - 1.0) +
                                p11 * (x * y - 2.0) + p02 * (y * y - 3.0));

    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }

  SECTION("saxpy-like") {
    auto x = zisa::XYZ{-3.4, 2.138, 0.0};

    auto pq = zisa::Poly2D<5, n_vars>(0.2 * p + q - 0.4 * p);

    auto exact = zisa::Cartesian<n_vars>{0.2 * p(x) + q(x) - 0.4 * p(x)};
    auto approx = pq(x);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }
}
