#include <zisa/testing/testing_framework.hpp>

#include <numeric>
#include <utility>
#include <vector>

#include <zisa/math/edge_rule.hpp>

TEST_CASE("Edge rule; basics") {

  // Poly *degree* -> n_points.
  std::vector<std::pair<zisa::int_t, zisa::int_t>> deg_n_points{
      {0, 1},{1, 1}, {2, 2}, {3, 2}, {4, 3}, {5, 3}, {6, 4}, {7, 4}};

  for (auto &&[deg, n_points] : deg_n_points) {

    auto qr = zisa::cached_edge_quadrature_rule(deg);

    REQUIRE(qr.points.size() == n_points);
    REQUIRE(qr.weights.size() == n_points);

    double total_weight
        = std::accumulate(qr.weights.begin(), qr.weights.end(), 0.0);

    REQUIRE(zisa::almost_equal(total_weight, 1.0, 1e-12));

    for (zisa::int_t k = 0; k < n_points; ++k) {
      REQUIRE(zisa::almost_equal(
          qr.points[k], -qr.points[n_points - 1 - k], 1e-12));

      REQUIRE(zisa::almost_equal(
          qr.weights[k], qr.weights[n_points - 1 - k], 1e-12));
    }
  }
}
