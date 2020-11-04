#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/few_points_cache.hpp>

TEST_CASE("FewPointsCache; API", "[math][cache]") {
  auto points = zisa::array<zisa::XYZ, 1>(10);
  points[0] = {0.52863914, 0.61846075, 0.00379484};
  points[1] = {0.96970712, 0.05154537, 0.01580046};
  points[2] = {0.51599433, 0.60140505, 0.20233078};
  points[3] = {0.11376839, 0.95724425, 0.37710804};
  points[4] = {0.85471734, 0.35608016, 0.72502333};
  points[5] = {0.69788378, 0.52131756, 0.84747612};
  points[6] = {0.31908423, 0.88898268, 0.80377053};
  points[7] = {0.60861228, 0.90587042, 0.17127816};
  points[8] = {0.24907181, 0.74370425, 0.5859984};
  points[9] = {0.68935398, 0.47523194, 0.26126917};

  auto cache = zisa::FewPointsCache<double>(points);

  auto f = [](const zisa::XYZ &x) { return zisa::norm(x); };

  cache.update(f);

  for (auto p : points) {
    REQUIRE(cache.get(p) == f(p));

    double eps = 1e-10;
    auto p_prime = zisa::XYZ{p[0] + eps, p[1] - eps, p[2] - eps};

    REQUIRE(cache.get(p) == f(p));
  }
}
