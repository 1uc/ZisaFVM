#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("GridVariables") {
  zisa::int_t n_cells = 10;
  zisa::int_t n_vars = 5;

  auto u = zisa::GridVariables(zisa::shape_t<2>{n_cells, n_vars});
  std::fill(u.begin(), u.end(), 0.0);

  auto ua = zisa::euler_var_t{1.0, 2.0, 3.0, 4.0, 5.0};

  u(2) = ua;
  for (zisa::int_t k = 0; k < n_vars; ++k) {
    REQUIRE(u(2, k) == ua(k));
  }

  ua = u(3);
  for (zisa::int_t k = 0; k < n_vars; ++k) {
    REQUIRE(ua(k) == 0.0);
  }
}
