#include <zisa/cuda/magma_batched.hpp>

namespace zisa {
namespace cuda {
namespace magma {

namespace internal {
ANY_DEVICE_INLINE
static double *extract_pointer(const array_const_view<double, 2, column_major> &a, int_t k) {
  return const_cast<double *>(&a(0, k));
}

ANY_DEVICE_INLINE
static double *extract_pointer(const array_const_view<double, 3, column_major> &a, int_t k) {
  return const_cast<double *>(&a(0, 0, k));
}

template <int NDIMS>
__global__ void fill_ptrs_kernel(array_view<double *, 1> a_ptrs,
                                 array_const_view<double, NDIMS, column_major> a) {
  auto tid = threadIdx.x;
  int_t i = tid;
  auto n_mat = a_ptrs.shape(0);

  while (i < n_mat) {
    a_ptrs[i] = extract_pointer(a, i);

    i += blockDim.x;
  }
}

template <int NDIMS>
void fill_ptrs_impl(const array_view<double *, 1> &a_ptrs,
                    const array_const_view<double, NDIMS, column_major> &a) {

  // TODO parallelize
  auto n_mat = a_ptrs.shape(0);
  for (int_t i = 0; i < n_mat; ++i) {
    a_ptrs[i] = extract_pointer(a, i);
  }
}

template <int NDIMS>
void fill_ptrs(const array_view<double *, 1> &a_ptrs,
               const array_const_view<double, NDIMS, column_major> &a) {

  auto mem_loc = zisa::memory_location(a_ptrs);
  assert(zisa::memory_location(a) == mem_loc);
  if (mem_loc == device_type::cpu) {
    internal::fill_ptrs_impl(a_ptrs, a);
  } else if (mem_loc == device_type::cuda) {
    internal::fill_ptrs_kernel<<<1024, 1>>>(a_ptrs, a);
    cudaDeviceSynchronize();
    ZISA_CHECK_CUDA;
  } else {
    LOG_ERR("Unknown memory location.");
  }
}
}

#define ZISA_FILL_PTRS_DEFN(NDIMS) \
void fill_ptrs(const array_view<double *, 1> &a_ptrs, \
const array_const_view<double, NDIMS, column_major> &a) { \
  internal::fill_ptrs(a_ptrs, a);  \
}

ZISA_FILL_PTRS_DEFN(2);
ZISA_FILL_PTRS_DEFN(3);

#undef ZISA_FILL_PTRS_DECL

}}}