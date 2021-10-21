#include <zisa/testing/testing_framework.hpp>

#include <zisa/cuda/magma_batched.hpp>
#include <zisa/memory/array.hpp>

TEST_CASE("MAGMA; extract pointers (cpu)", "[magma]") {
  zisa::int_t n_mat = 6;
  zisa::int_t m = 4;
  zisa::int_t n = 5;
  auto a = zisa::array<double, 3, zisa::column_major>(zisa::shape_t<3>{m, n, n_mat});
  auto a_ptrs = zisa::array<double *, 1>(n_mat);

  zisa::cuda::magma::fill_ptrs(a_ptrs, a);

  for(zisa::int_t k = 0; k < n_mat; ++k) {
    REQUIRE(a_ptrs(k) == a.raw() + k*m*n);
  }
}

TEST_CASE("MAGMA; extract pointers (cuda)", "[magma]") {
  zisa::int_t n_mat = 6;
  zisa::int_t m = 4;
  zisa::int_t n = 5;
  auto a = zisa::array<double, 3, zisa::column_major>(zisa::shape_t<3>{m, n, n_mat}, zisa::device_type::cuda);
  auto a_ptrs_device = zisa::array<double *, 1>(n_mat, zisa::device_type::cuda);

  zisa::cuda::magma::fill_ptrs(a_ptrs_device, a);

  auto a_ptrs = zisa::array<double *, 1>(n_mat, zisa::device_type::cpu);
  zisa::copy(a_ptrs, a_ptrs_device);

  for(zisa::int_t k = 0; k < n_mat; ++k) {
    REQUIRE(a_ptrs(k) == a.raw() + k*m*n);
  }
}
