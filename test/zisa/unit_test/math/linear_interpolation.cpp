
#include <zisa/math/comparison.hpp>
#include <zisa/math/linear_interpolation.hpp>
#include <zisa/testing/testing_framework.hpp>

namespace zisa {

TEST_CASE("Linear interpolataion", "[math]") {
  auto points = array<double, 1>(3ul);
  auto values = array<double, 1>(3ul);

  points[0] = 0.0;
  points[1] = 1.0;
  points[2] = 3.0;

  values[0] = 1.0;
  values[1] = 2.0;
  values[2] = 6.0;

  auto interp = zisa::NonUniformLinearInterpolation(points, values);

  double atol = 1e-13;

  REQUIRE(zisa::almost_equal(interp(0.0), 1.0, atol));
  REQUIRE(zisa::almost_equal(interp(0.5), 1.5, atol));
  REQUIRE(zisa::almost_equal(interp(1.0), 2.0, atol));
  REQUIRE(zisa::almost_equal(interp(2.0), 4.0, atol));
  REQUIRE(zisa::almost_equal(interp(2.5), 5.0, atol));
}

}