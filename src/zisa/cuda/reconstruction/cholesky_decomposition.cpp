#include <zisa/cuda/reconstruction/cholesky_decomposition.hpp>

namespace zisa::cuda::magma {

CholeskyDecomposition::CholeskyDecomposition(
    std::shared_ptr<MAGMAQueue> magma_context_,
    zisa::array<double, 3, zisa::column_major> A_batched_)
    : magma_queue(std::move(magma_context_)) {
  allocate(std::move(A_batched_));
  decompose();
}
void CholeskyDecomposition::solve(array<double, 3, column_major> &rhs_batched) {
  zisa::cuda::magma::fill_ptrs(rhs_ptrs, rhs_batched);
  solve_cholesky_inplace(
      ll_batched, ll_ptrs, rhs_batched, rhs_ptrs, *magma_queue);
}

void CholeskyDecomposition::allocate(array<double, 3, column_major> A_batched) {
  auto mem_loc = zisa::memory_location(A_batched);
  int_t n_mats = A_batched.shape(2);

  ll_batched = std::move(A_batched);
  ll_ptrs = array<double *, 1>(n_mats, mem_loc);
  fill_ptrs(ll_ptrs, ll_batched);

  rhs_ptrs = array<double *, 1>(n_mats, mem_loc);
  info_batched = array<magma_int_t, 1>(n_mats);
}

void CholeskyDecomposition::decompose() {
  decompose_cholesky_inplace(ll_batched, ll_ptrs, info_batched, *magma_queue);
}

void decompose_cholesky_inplace(
    const array_view<double, 3, column_major> &ll_batched,
    const array_view<double *, 1> &ll_ptrs,
    const array_view<magma_int_t, 1> &info_batched,
    const MAGMAQueue &magma_queue) {

  assert(zisa::memory_location(ll_batched) == device_type::cuda);
  assert(zisa::memory_location(ll_ptrs) == device_type::cuda);

  int_t n_rows = ll_batched.shape(0);
  int_t n_cols = ll_batched.shape(1);
  int_t n_mats = ll_batched.shape(2);

  assert(n_rows == n_cols);

  auto uplo = MagmaLower;
  auto n = integer_cast<magma_int_t>(n_cols);
  auto dA_array = ll_ptrs.raw();
  auto ldda = integer_cast<magma_int_t>(n_rows);
  auto info_array = info_batched.raw();
  auto batchCount = integer_cast<magma_int_t>(n_mats);
  auto queue = magma_queue.queue();

  auto status = magma_dpotrf_batched(
      uplo, n, dA_array, ldda, info_array, batchCount, queue);
  magma_queue.sync();

  LOG_ERR_IF(status != MAGMA_SUCCESS, "An error in MAGMA occurred.");
}

void solve_cholesky_inplace(
    const array_view<double, 3, column_major> &ll_batched,
    const array_view<double *, 1> &ll_ptrs,
    const array_view<double, 3, column_major> &rhs_batched,
    const array_view<double *, 1> &rhs_ptrs,
    const MAGMAQueue &magma_queue) {

  int_t n_rows = ll_batched.shape(0);
  int_t n_mats = ll_batched.shape(2);
  int_t n_rhs = rhs_batched.shape(1);

  auto side = MagmaLeft;
  auto uplo = MagmaLower;
  auto diag = MagmaNonUnit;
  auto m = integer_cast<magma_int_t>(n_rows);
  auto n = integer_cast<magma_int_t>(n_rhs);
  auto alpha = 1.0;
  auto dA_array = ll_ptrs.raw();
  auto ldda = m;

  auto dB_array = rhs_ptrs.raw();
  auto lddb = m;
  auto batchCount = integer_cast<magma_int_t>(n_mats);
  auto queue = magma_queue.queue();

  // Solve:
  // L * C = B

  // clang-format off
  magmablas_dtrsm_batched(
    side, uplo, MagmaNoTrans, diag,
    m, n, alpha,
    dA_array, ldda,
    dB_array, lddb,
    batchCount, queue
  );
  // clang-format on

  // Solve:
  // L^T * X = C

  // clang-format off
  magmablas_dtrsm_batched(
    side, uplo, MagmaTrans, diag,
    m, n, alpha,
    dA_array, ldda,
    dB_array, lddb,
    batchCount, queue
  );
  // clang-format on
}

}