#include <catch/catch.hpp>

#include "zisa/memory/array.hpp"
#include "zisa/memory/array_view.hpp"
#include "zisa/memory/column_major.hpp"

using namespace zisa;

bool implicit_conversion(
    const array_const_view<double, 3, column_major> &view) {

  return view.raw() != nullptr;
}

TEST_CASE("array; basics") {
  auto a = array<double, 3>({3ul, 3ul, 2ul}, device_type::cpu);
  auto b = array<double, 3>({3ul, 3ul, 2ul}, device_type::cpu);

  REQUIRE(a.raw() != static_cast<const decltype(b) &>(b).raw());
  REQUIRE(a.shape() == b.shape());

  auto b_view = array_view<double, 3>(b);
  auto bb_view = array_view<double, 3>(shape(b), raw_ptr(b));
  auto bc_view = array_const_view<double, 3>(shape(b), raw_ptr(b));

  REQUIRE(b.raw() == b_view.raw());

  b(1, 0, 1) = -42.0;
  REQUIRE(b_view(1, 0, 1) == -42.0);

  b_view(0, 0, 0) = 42.0;
  REQUIRE(b(0, 0, 0) == 42.0);

  REQUIRE(implicit_conversion(a));
}
