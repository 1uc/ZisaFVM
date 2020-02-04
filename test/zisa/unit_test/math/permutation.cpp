#include <zisa/testing/testing_framework.hpp>

#include <zisa/math/permutation.hpp>


TEST_CASE("Permuation", "[math]") {

  // The example is:
  //  (0) (1 2 4 5) (3 6)
  // which is
  //  0 1 2 3 4 5 6
  //  0 2 4 6 5 1 3

  std::vector<zisa::int_t> sigma{0, 2, 4, 6, 5, 1, 3};

  auto permutation = zisa::factor_permutation(sigma);
  const auto &cycles = permutation.cycles;

  CHECK(cycles.size() == 3);

  CHECK(cycles(0).size() == 1);
  CHECK(cycles(0)(0) == 0);

  CHECK(cycles(1).size() == 4);
  CHECK(cycles(1)(0) == 1);
  CHECK(cycles(1)(1) == 2);
  CHECK(cycles(1)(2) == 4);
  CHECK(cycles(1)(3) == 5);

  CHECK(cycles(2).size() == 2);
  CHECK(cycles(2)(0) == 3);
  CHECK(cycles(2)(1) == 6);

  std::vector<double> original{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
  std::vector<double> to_permute = original;
  std::vector<double> permutated{0.0, 20.0, 40.0, 60.0, 50.0, 10.0, 30.0};

  zisa::apply_permutation(zisa::array_view(to_permute), permutation);
  CHECK(to_permute == permutated);

  zisa::reverse_permutation(zisa::array_view(to_permute), permutation);
  CHECK(to_permute == original);
}