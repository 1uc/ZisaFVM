#ifndef MAGMA_BATCHED_HPP
#define MAGMA_BATCHED_HPP

#include <magma_v2.h>

#include <zisa/memory/array_view.hpp>

namespace zisa {
namespace cuda {
namespace magma {

#define ZISA_FILL_PTRS_DECL(NDIMS) \
void fill_ptrs(const array_view<double *, 1> &a_ptrs, \
               const array_const_view<double, NDIMS, column_major> &a);

ZISA_FILL_PTRS_DECL(2)
ZISA_FILL_PTRS_DECL(3)

#undef ZISA_FILL_PTRS_DECL

}
}
}
#endif // MAGMA_BATCHED_HPP