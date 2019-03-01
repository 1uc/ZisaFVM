#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/barycentric.hpp>

TEST_CASE("Barycentric; basic API", "[math]") {
  auto tri = zisa::Triangle{{0.461024, 0.432542, 0.0},
                            {0.558699, 0.333827, 0.0},
                            {0.554918, 0.426801, 0.0}};
  auto x_center = zisa::XYZ{{0.52488, 0.397723, 0.0}};

  auto lambda = zisa::Barycentric(tri, x_center);

  for (int k = 0; k < 3; ++k) {
    INFO(string_format("lambda[%d] = %e", k, lambda[k]));
    REQUIRE(zisa::almost_equal(lambda[k], 1.0 / 3.0, 1e-4));
  }
}
