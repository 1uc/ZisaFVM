// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/poly2d.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/math/triangle.hpp>

TEST_CASE("poly_dof-2d", "[math][poly]") {
  REQUIRE(zisa::poly_dof<2>(0) == 1);
  REQUIRE(zisa::poly_dof<2>(1) == 3);
  REQUIRE(zisa::poly_dof<2>(2) == 6);
  REQUIRE(zisa::poly_dof<2>(3) == 10);
  REQUIRE(zisa::poly_dof<2>(4) == 15);
}

TEST_CASE("poly_dof-3d", "[math][poly]") {
  REQUIRE(zisa::poly_dof<3>(0) == 1);
  REQUIRE(zisa::poly_dof<3>(1) == 4);
  REQUIRE(zisa::poly_dof<3>(2) == 10);
  REQUIRE(zisa::poly_dof<3>(3) == 20);
  REQUIRE(zisa::poly_dof<3>(4) == 35);
}

TEST_CASE("PolyIndexRange<2>", "[math][poly]") {
  zisa::int_t max_degree = 3;
  auto index_range = zisa::PolyIndexRange<2>(max_degree);

  zisa::int_t degree = 0;
  zisa::int_t count = 0;
  for (auto [i, j] : index_range) {
    if (i + j > degree) {
      ++degree;
    }

    REQUIRE(i + j == degree);
    REQUIRE(zisa::poly_index(i, j) == count);
    ++count;
  }

  REQUIRE(degree == max_degree);
  REQUIRE(count == zisa::poly_dof<2>(max_degree));
}

TEST_CASE("PolyIndexRange<3>", "[math][poly]") {
  zisa::int_t max_degree = 3;
  auto index_range = zisa::PolyIndexRange<3>(max_degree);

  zisa::int_t degree = 0;
  zisa::int_t count = 0;
  for (auto [i, j, k] : index_range) {
    if (i + j + k > degree) {
      ++degree;
    }

    REQUIRE(i + j + k == degree);
    REQUIRE(zisa::poly_index(i, j, k) == count);
    ++count;
  }

  REQUIRE(max_degree == degree);
  REQUIRE(count == zisa::poly_dof<3>(max_degree));
}

template <int NDIMS>
void check_poly_degree_right_inverse() {
  for (zisa::int_t n_coeffs = 1; n_coeffs < 25; ++n_coeffs) {
    auto deg = zisa::poly_degree<NDIMS>(n_coeffs);

    INFO(string_format("n_coeffs = %d \n", n_coeffs));
    REQUIRE(zisa::poly_dof<NDIMS>(deg) <= n_coeffs);
    REQUIRE(zisa::poly_dof<NDIMS>(deg + 1) > n_coeffs);
  }
}

template <int NDIMS>
void check_poly_degree_left_inverse() {
  for (zisa::int_t deg = 0; deg < 5; ++deg) {
    INFO(string_format("deg = %d \n", deg));
    REQUIRE(zisa::poly_degree<NDIMS>(zisa::poly_dof<NDIMS>(deg)) == deg);
  }
}

TEST_CASE("poly_degree", "[math][poly]") {

  SECTION("right-inverse-2d") { check_poly_degree_right_inverse<2>(); }
  SECTION("right-inverse-3d") { check_poly_degree_right_inverse<3>(); }

  SECTION("left-inverse-2d") { check_poly_degree_left_inverse<2>(); }
  SECTION("left-inverse-3d") { check_poly_degree_left_inverse<3>(); }
}

TEST_CASE("Poly2D; examples", "[math][poly]") {
  int n_dims = 2;
  using Poly = zisa::PolyND<zisa::poly_dof<2>(5), 1>;
  auto p = Poly({1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                {0.0, 0.0, 0.0, 1.0, 2.0, 3.0},
                zisa::XYZ::zeros(),
                1.0,
                n_dims);

  auto q
      = Poly({1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}, zisa::XYZ::zeros(), 1.0, n_dims);

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
    auto exact = 1.0 + 2.0 * x + 3.0 * y + 4.0 * (x * x - 1.0)
                 + 5.0 * (x * y - 2.0) + 6.0 * (y * y - 3.0);

    INFO(string_format("%e != %e [%e]\n", approx[0], exact, approx - exact));
    REQUIRE(zisa::almost_equal(approx[0], exact, 1e-14));
  }

  SECTION("saxpy-like") {
    auto x = zisa::XYZ{-3.4, 2.138, 0.0};

    auto pq = Poly(0.2 * p + q - 0.4 * p);

    auto exact = zisa::Cartesian<1>{0.2 * p(x) + q(x) - 0.4 * p(x)};
    auto approx = pq(x);

    INFO(string_format("%e != %e [%e]\n", approx[0], exact[0], approx - exact));
    REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
  }

  SECTION("assignment") {
    auto x = zisa::XYZ{-3.4, 2.138, 0.0};

    auto tmp = Poly(p);

    SECTION("+=") {
      tmp += q;

      auto exact = zisa::Cartesian<1>{p(x) + q(x)};
      auto approx = tmp(x);

      INFO(string_format(
          "%e != %e [%e]\n", approx[0], exact[0], approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }

    SECTION("+=") {
      tmp -= q;

      auto exact = zisa::Cartesian<1>{p(x) - q(x)};
      auto approx = zisa::Cartesian<1>{tmp(x)};

      INFO(string_format(
          "%e != %e [%e]\n", approx[0], exact[0], approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }

    SECTION("*=") {
      tmp *= 2.0;

      auto exact = zisa::Cartesian<1>{2.0 * p(x)};
      auto approx = zisa::Cartesian<1>{tmp(x)};

      INFO(string_format(
          "%e != %e [%e]\n", approx[0], exact[0], approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }

    SECTION("/=") {
      tmp /= 2.0;

      auto exact = zisa::Cartesian<1>{0.5 * p(x)};
      auto approx = zisa::Cartesian<1>{tmp(x)};

      INFO(string_format(
          "%e != %e [%e]\n", approx[0], exact[0], approx - exact));
      REQUIRE(zisa::almost_equal(approx[0], exact[0], 1e-14));
    }
  }

  SECTION("quadrature") {
    auto tri = zisa::reference_triangle();

    auto pbar = zisa::quadrature(p, tri, 4);
    REQUIRE(pbar[0] < 0.0);
  }
}

TEST_CASE("Poly3D<., 2>; examples", "[math][poly]") {
  int n_dims = 3;
  constexpr int n_vars = 2;
  using Poly = zisa::PolyND<zisa::poly_dof<3>(5), n_vars>;

  auto p = Poly(2,
                {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                zisa::XYZ::zeros(),
                1.0,
                n_dims);

  auto p000 = zisa::Cartesian<n_vars>{1.0, -1.0};

  auto p100 = zisa::Cartesian<n_vars>{2.0, -1.0};
  auto p010 = zisa::Cartesian<n_vars>{-1.0, -1.0};
  auto p001 = zisa::Cartesian<n_vars>{-3.0, -1.0};

  auto p200 = zisa::Cartesian<n_vars>{2.0, -3.0};
  auto p110 = zisa::Cartesian<n_vars>{4.0, -2.0};
  auto p101 = zisa::Cartesian<n_vars>{2.0, -5.0};
  auto p020 = zisa::Cartesian<n_vars>{1.0, -3.0};
  auto p011 = zisa::Cartesian<n_vars>{1.0, -3.0};
  auto p002 = zisa::Cartesian<n_vars>{3.0, -1.0};

  // clang-format off
  auto p_coeffs = std::vector<zisa::Cartesian<n_vars>>{
         p000,
         p100, p010, p001,
         p200, p110, p101, p020, p011, p002
  };
  // clang-format on

  std::copy((double *)&p_coeffs[0],
            (double *)&p_coeffs[0] + p.n_vars() * p.dof(),
            p.coeffs_ptr());

  auto q = Poly(1, {0.0, 0.0, 0.0, 0.0}, zisa::XYZ::zeros(), 1.0, n_dims);

  auto q000 = zisa::Cartesian<n_vars>{1.0, -1.0};
  auto q100 = zisa::Cartesian<n_vars>{-1.0, -1.0};
  auto q010 = zisa::Cartesian<n_vars>{2.0, -1.0};
  auto q001 = zisa::Cartesian<n_vars>{3.0, -1.0};

  auto q_coeffs = std::vector<zisa::Cartesian<n_vars>>{q000, q100, q010, q001};
  std::copy((double *)&q_coeffs[0],
            (double *)&q_coeffs[0] + q.n_vars() * q.dof(),
            q.coeffs_ptr());

  SECTION("degree") { REQUIRE(p.degree() == 2); }
  SECTION("max_degree") { REQUIRE(p.max_degree() == 5); }
  SECTION("n_vars") { REQUIRE(p.n_vars() == n_vars); }

  SECTION("p{} == 0") {
    auto p = std::make_shared<Poly>();
    REQUIRE((*p)(zisa::XYZ{2.0, -1.0, 0.0}) == zisa::Cartesian<n_vars>(0.0));
  }

  SECTION("p(x)") {
    double x = -3.0, y = 2.0, z = 0.2;

    auto approx = p(zisa::XYZ{x, y, z});

    // clang-format off
    auto exact
        = zisa::Cartesian<n_vars>(
              p000
            + p100 * x + p010 * y + p001 * z
            + p200 * (x*x - 1.0) + p110 * (x*y - 2.0) + p101 * (x*z - 3.0)
            + p020 * (y*y - 4.0) + p011 * (y*z - 5.0) + p002 * (z*z - 6.0));
    // clang-format on

    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }

  SECTION("saxpy-like") {
    auto x = zisa::XYZ{-3.4, 2.138, 0.2};

    auto pq = Poly(0.2 * p + q - 0.4 * p);

    auto exact = zisa::Cartesian<n_vars>{0.2 * p(x) + q(x) - 0.4 * p(x)};
    auto approx = pq(x);

    REQUIRE(zisa::almost_equal(approx, exact, 1e-14));
  }
}
