#include <catch/catch.hpp>

#include <algorithm>
#include <iostream>
#include <memory>

#include "zisa/memory/contiguous_memory.hpp"
#include "zisa/memory/host_contiguous_memory.hpp"

template <class T>
zisa::contiguous_memory<bool> check_equal(const zisa::contiguous_memory<T> &a,
                                          const zisa::contiguous_memory<T> &b) {
  assert(a.size() == b.size());

  auto is_equal = zisa::contiguous_memory<bool>(a.size());

  for (int i = 0; i < a.size(); ++i) {
    is_equal[i] = (a[i] == b[i]);
  }

  return is_equal;
}

void check_copyless_const_ref(const zisa::contiguous_memory<double> &a,
                              void *ptr) {
  REQUIRE(a.raw() == ptr);
}

void check_copyless_ref(zisa::contiguous_memory<double> &a, void *ptr) {
  REQUIRE(a.raw() == ptr);
}

template <class SRC>
void check_copy_construction(const SRC &src) {
  auto cpy = src;

  REQUIRE(cpy.raw() != src.raw());

  auto is_equal = check_equal(cpy, src);
  auto is_good
      = std::all_of(is_equal.begin(), is_equal.end(), [](bool x) { return x; });

  REQUIRE(is_good);
}

template <class SRC>
void check_move_construction(SRC src) {
  auto ptr = src.raw();

  auto cpy = std::move(src);

  REQUIRE(cpy.raw() == ptr);
  REQUIRE(src.raw() == nullptr);
}

TEST_CASE("contiguous_memory") {
  int n_elements = 15;

  auto a = zisa::contiguous_memory<double>(n_elements);
  std::fill(a.begin(), a.end(), 0.0);

  SECTION("copy-less upcast -- ref") {
    SECTION("host") { check_copyless_ref(a, a.raw()); }
  }

  SECTION("copy-less upcast -- const ref") {
    SECTION("host") { check_copyless_const_ref(a, a.raw()); }
  }

  SECTION("copy construction") {
    SECTION("host") { check_copy_construction(a); }
  }

  SECTION("move construction") {
    SECTION("host") { check_move_construction(a); }
  }
}

TEST_CASE("contiguous_memory; initialization", "[memory]") {
  zisa::int_t n_elements = 15;
  zisa::int_t n_arrays = 5;

  auto elem = zisa::contiguous_memory<double>(n_elements);
  std::fill(elem.begin(), elem.end(), 42.0);

  auto seq = zisa::contiguous_memory<zisa::contiguous_memory<double>>(n_arrays);

  for (zisa::int_t i = 0; i < n_arrays; ++i) {
    REQUIRE(seq[i].raw() == nullptr);
  }

  for (zisa::int_t i = 0; i < n_arrays; ++i) {
    elem[0] = 42.0 + 0.1 * i;
    seq[i] = elem;
  }

  for (zisa::int_t i = 0; i < n_arrays; ++i) {
    elem[0] = 42.0 + 0.1 * i;
    REQUIRE(std::equal(seq[i].begin(), seq[i].end(), elem.begin(), elem.end()));
  }
}
