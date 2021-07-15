// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <vector>

#include <zisa/math/symmetric_choices.hpp>
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("Symmetric choices", "[math]") {

  std::vector<zisa::int_t> vi{3, 22, 1, 6};

  auto choices = zisa::strict_symmetric_choices(vi);

  auto it = choices.begin();
  auto end_it = choices.end();

  auto check = [&end_it, &it](zisa::int_t i, zisa::int_t j) {
    auto [k, l] = *it;
    REQUIRE(i == k);
    REQUIRE(j == l);
    REQUIRE(it != end_it);

    ++it;
  };

  check(3, 22);
  check(3, 1);
  check(3, 6);

  check(22, 1);
  check(22, 6);

  check(1, 6);

  REQUIRE(it == end_it);
}

TEST_CASE("Symmetric choices; empty", "[math][runme]") {

  std::vector<zisa::int_t> vi{};

  auto choices = zisa::strict_symmetric_choices(vi);

  auto it = choices.begin();
  auto end_it = choices.end();

  REQUIRE(it == end_it);
}
