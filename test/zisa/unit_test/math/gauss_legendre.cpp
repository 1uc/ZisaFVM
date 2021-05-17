#include <zisa/testing/testing_framework.hpp>

#include "zisa/math/gauss_legendre.hpp"

namespace zisa {

TEST_CASE("GaussLegendre; five_point", "[math]") {
  int n_qp = 5;
  double points[5]
      = {-1.0 / 3.0 * zisa::sqrt(5.0 + 2.0 * zisa::sqrt(10.0 / 7.0)),
         -1.0 / 3.0 * zisa::sqrt(5.0 - 2.0 * zisa::sqrt(10.0 / 7.0)),
         0.0,
         1.0 / 3.0 * zisa::sqrt(5.0 - 2.0 * zisa::sqrt(10.0 / 7.0)),
         1.0 / 3.0 * zisa::sqrt(5.0 + 2.0 * zisa::sqrt(10.0 / 7.0))};

  double weights[5] = {(322.0 - 13.0 * zisa::sqrt(70.0)) / 900.0,
                       (322.0 + 13.0 * zisa::sqrt(70.0)) / 900.0,
                       128.0 / 225.0,
                       (322.0 + 13.0 * zisa::sqrt(70.0)) / 900.0,
                       (322.0 - 13.0 * zisa::sqrt(70.0)) / 900.0};

  GaussLegendre<5> gl;

  for (int i = 0; i < n_qp; ++i) {
    INFO(string_format("Failed on i = %d.", i));
    REQUIRE(almost_equal(gl.points[i], points[i], 2e-14));

    INFO(string_format("Failed on i = %d.", i));
    REQUIRE(almost_equal(gl.weights[i], weights[i], 2e-14));
  }
}

TEST_CASE("GaussLegendre; four_point", "[math]") {
  constexpr int n_qp = 4;
  double points[n_qp]
      = {-zisa::sqrt(3.0 / 7.0 + 2.0 / 7.0 * zisa::sqrt(6.0 / 5.0)),
         -zisa::sqrt(3.0 / 7.0 - 2.0 / 7.0 * zisa::sqrt(6.0 / 5.0)),
         zisa::sqrt(3.0 / 7.0 - 2.0 / 7.0 * zisa::sqrt(6.0 / 5.0)),
         zisa::sqrt(3.0 / 7.0 + 2.0 / 7.0 * zisa::sqrt(6.0 / 5.0))};

  double weights[n_qp] = {(18.0 - zisa::sqrt(30)) / 36.0,
                          (18.0 + zisa::sqrt(30)) / 36.0,
                          (18.0 + zisa::sqrt(30)) / 36.0,
                          (18.0 - zisa::sqrt(30)) / 36.0};

  GaussLegendre<n_qp> gl;

  for (int i = 0; i < n_qp; ++i) {
    INFO(string_format("Failed on i = %d.", i));
    REQUIRE(almost_equal(gl.points[i], points[i], 2e-14));

    INFO(string_format("Failed on i = %d.", i));
    REQUIRE(almost_equal(gl.weights[i], weights[i], 2e-14));
  }
}

TEST_CASE("GaussLegendre; check_roots_five", "[math]") {
  constexpr int n_qp = 5;
  double points[3]
      = {0.0,
         1.0 / 3.0 * zisa::sqrt(5.0 - 2.0 * zisa::sqrt(10.0 / 7.0)),
         1.0 / 3.0 * zisa::sqrt(5.0 + 2.0 * zisa::sqrt(10.0 / 7.0))};

  auto fn = FourierNewton();
  auto p = fn.fourier_poly<n_qp>();

  REQUIRE(almost_equal(p.at(n_qp),
                       zisa::sqrt(5.5) * double(3 * 5 * 7 * 9)
                           / (2 * 2 * 2 * 2 * 5.0 * 4 * 3 * 2),
                       1e-12));

  for (int i = 0; i < 3; ++i) {
    INFO(string_format("Failed on i = %d.", i));
    REQUIRE(almost_equal(p(::acos(zisa::abs(points[i]))), 0.0, 1e-8));
  }
}

TEST_CASE("GaussLegendre; check_roots_four", "[math]") {
  constexpr int n_qp = 4;
  double points[2]
      = {zisa::sqrt(3.0 / 7.0 - 2.0 / 7.0 * zisa::sqrt(6.0 / 5.0)),
         zisa::sqrt(3.0 / 7.0 + 2.0 / 7.0 * zisa::sqrt(6.0 / 5.0))};

  auto fn = FourierNewton();
  auto p = fn.fourier_poly<n_qp>();

  double a_nn
      = zisa::sqrt(4.5) * double(3 * 5 * 7) / double(2 * 2 * 2 * 4 * 3 * 2);
  REQUIRE(almost_equal(p.at(4), a_nn, 1e-12));
  REQUIRE(almost_equal(p.at(2), 4.0 / (8.0 - 1) * a_nn, 1e-12));
  REQUIRE(
      almost_equal(p.at(0),
                   0.5 * 3.0 / 2.0 * 4.0 * 3.0 / ((8.0 - 1) * (8.0 - 3)) * a_nn,
                   1e-12));

  for (int i = 0; i < 2; ++i) {
    INFO(string_format("Failed on i = %d.", i));
    REQUIRE(almost_equal(p(::acos(zisa::abs(points[i]))), 0.0, 1e-8));
  }
}

} // namespace zisa
