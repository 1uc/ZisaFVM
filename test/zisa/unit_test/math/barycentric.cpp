// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/barycentric.hpp>

TEST_CASE("Barycentric2D; is_inside", "[math][barycentric]") {

  auto tri = zisa::reference_triangle();

  REQUIRE(zisa::is_inside(zisa::Barycentric2D(tri, zisa::XYZ{0.1, 0.1, 0.0})));
  REQUIRE(
      !zisa::is_inside(zisa::Barycentric2D(tri, zisa::XYZ{-0.1, 0.1, 0.0})));
  REQUIRE(!zisa::is_inside(zisa::Barycentric2D(tri, zisa::XYZ{1.1, 1.1, 0.0})));
}

TEST_CASE("Barycentric2D; basic API", "[math][barycentric]") {
  auto tri = zisa::Triangle{{0.461024, 0.432542, 0.0},
                            {0.558699, 0.333827, 0.0},
                            {0.554918, 0.426801, 0.0}};

  auto x = zisa::XYZ{{0.588, 0.897723, 0.0}};
  auto lambda = zisa::Barycentric2D(tri, x);
  auto xx = zisa::coord(tri, lambda);

  for (int k = 0; k < 3; ++k) {
    INFO(string_format("lambda[%d] = %e", k, lambda[k]));
    REQUIRE(zisa::almost_equal(x[k], xx[k], 1e-10));
  }
}

TEST_CASE("Barycentric3D; basic API", "[math][barycentric]") {
  auto tet = zisa::Tetrahedron{
      {0.554918, 0.426801, 0.1},
      {0.461024, 0.432542, 0.05},
      {0.558699, 0.333827, 0.0012},
      {0.504918, 0.406801, 1.0},
  };
  auto x = zisa::XYZ{{0.50088, 0.3923, 0.2}};
  auto lambda = zisa::Barycentric3D(tet, x);
  auto xx = zisa::coord(tet, lambda);

  for (int k = 0; k < 3; ++k) {
    INFO(string_format("lambda = %e, x = %e, xx = %e, dx = %e, k =  %d",
                       lambda[k],
                       x[k],
                       xx[k],
                       x[k] - xx[k],
                       k));
    REQUIRE(zisa::almost_equal(x[k], xx[k], 1e-10));
  }
}
