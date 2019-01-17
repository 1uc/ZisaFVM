#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/triangle.hpp>
#include <zisa/math/basic_functions.hpp>

TEST_CASE("avg_moment", "[.]") {
  auto tri = zisa::Triangle({0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0});
  zisa::int_t deg = 3;

  auto exact = [&tri](int k, int l) {
    auto vol = tri.volume;
    auto m = zisa::Gamma(k + 2) * zisa::Gamma(l + 1)
             / ((k + 1) * zisa::Gamma(k + l + 3));

    return m / vol;
  };

  for (int k = 1; k < deg; ++k) {
    for (int l = 0; l < k; ++l) {
      auto approx = avg_moment(tri, k, l, 2);

      INFO(string_format("[%d, %d] %e  !=  %e \n", k, l, approx, exact(k, l)));
      REQUIRE(zisa::almost_equal(approx, exact(k, l), 0.2));
    }
  }
}
