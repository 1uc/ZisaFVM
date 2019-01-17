#include <iostream>
#include <memory>

#include "zisa/memory/contiguous_memory_base.hpp"
#include <zisa/testing/testing_framework.hpp>

TEST_CASE("array; STL allocator") {
  int n_elements = 15;

  auto alloc = std::allocator<double>();
  auto a = zisa::contiguous_memory_base<double, decltype(alloc)>(n_elements,
                                                                 alloc);
  auto b = zisa::contiguous_memory_base<double, decltype(alloc)>(n_elements,
                                                                 alloc);

  a[2] = 42;

  SECTION("move constructor") {
    auto ptr = a.raw();

    b = std::move(a);

    REQUIRE(b.raw() == ptr);
    REQUIRE(a.raw() == nullptr);
  }
}
