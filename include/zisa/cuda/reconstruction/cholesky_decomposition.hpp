#ifndef CHOLESKY_DECOMPOSITION_HPP
#define CHOLESKY_DECOMPOSITION_HPP

#include <zisa/config.hpp>
#include <zisa/cuda/magma_batched.hpp>
#include <zisa/cuda/magma_context.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/column_major.hpp>

namespace zisa {
namespace cuda {
namespace magma {

class CholeskyDecomposition {
public:
  /// Create a Cholesky based linear solver.
  /**
   * @param magma_context
   * @param A_batched has shape (n, n, k) representing 'k' square matrices n by
   *                  n matrices and A_batched(i, j, k) is the element A_{i,j}
   *                  in the k-th matrix.
   */
  CholeskyDecomposition(std::shared_ptr<MAGMAQueue> magma_context_,
                        array<double, 3, column_major> A_batched_);

  void solve(array<double, 3, column_major> &rhs);

private:
  void allocate(array<double, 3, column_major> A_batched);

  void decompose();

private:
  array<double, 3, column_major> ll_batched;
  array<double *, 1> ll_ptrs;

  array<double *, 1> rhs_ptrs;
  array<magma_int_t, 1> info_batched;

  std::shared_ptr<MAGMAQueue> magma_queue;
};

void solve_cholesky_inplace(
    const array_view<double, 3, column_major> &ll_batched,
    const array_view<double *, 1> &ll_ptrs,
    const array_view<double, 3, column_major> &rhs_batched,
    const array_view<double *, 1> &rhs_ptrs,
    const MAGMAQueue &magma_queue);

void decompose_cholesky_inplace(
    const array_view<double, 3, column_major> &ll_batched,
    const array_view<double *, 1> &ll_ptrs,
    const array_view<magma_int_t, 1> &info_batched,
    const MAGMAQueue &magma_queue);

}
}
}

#endif // CHOLESKY_DECOMPOSITION_HPP
