#include <zisa/testing/testing_framework.hpp>

#include <Eigen/Dense>
#include <zisa/cuda/reconstruction/cholesky_decomposition.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

static std::vector<Eigen::MatrixXd>
make_random_matrices(int_t n_rows, int_t n_cols, int_t n_mat) {
  auto matrices = std::vector<Eigen::MatrixXd>();
  matrices.reserve(n_mat);

  for (int_t k = 0; k < n_mat; ++k) {
    matrices.emplace_back(Eigen::MatrixXd::Random(n_rows, n_cols));
  }

  return matrices;
}

static std::vector<Eigen::MatrixXd> make_spd_matrices(int_t n_cols,
                                                      int_t n_mat) {
  auto matrices = make_random_matrices(n_cols, n_cols, n_mat);

  for (int_t k = 0; k < n_mat; ++k) {
    matrices[k] = matrices[k].transpose() * matrices[k];
  }

  return matrices;
}

array<double, 3, column_major>
convert_to_batched(const std::vector<Eigen::MatrixXd> &matrices) {

  auto n_mat = integer_cast<int_t>(matrices.size());
  auto n_cols = integer_cast<int_t>(matrices[0].cols());
  auto n_rows = integer_cast<int_t>(matrices[0].rows());

  auto A_batched
      = array<double, 3, column_major>(shape_t<3>{n_rows, n_cols, n_mat});
  for (int_t k = 0; k < n_mat; ++k) {
    for (int_t j = 0; j < n_cols; ++j) {
      for (int_t i = 0; i < n_rows; ++i) {
        A_batched(i, j, k) = matrices[k](i, j);
      }
    }
  }

  return A_batched;
}

void check_llt_against_eigen(int_t n_mat, int_t n_cols) {
  auto matrices = make_spd_matrices(n_cols, n_mat);

  auto llts = std::vector<Eigen::MatrixXd>{};
  llts.reserve(n_mat);

  for (int_t k = 0; k < n_mat; ++k) {
    llts.emplace_back(matrices[k].llt().matrixL());
  }

  auto ll_batched_host = convert_to_batched(matrices);
  auto ll_batched_device = array<double, 3, column_major>(
      ll_batched_host.shape(), device_type::cuda);
  zisa::copy(ll_batched_device, ll_batched_host);

  auto ll_ptrs_device = array<double *, 1>(n_mat, device_type::cuda);
  auto ll_ptrs_host = array<double *, 1>(n_mat, device_type::cpu);

  zisa::cuda::magma::fill_ptrs(ll_ptrs_device, ll_batched_device);
  auto info_batched = array<magma_int_t, 1>(n_mat, device_type::cuda);

  auto magma_queue = zisa::cuda::magma::make_default_queue();

  decompose_cholesky_inplace(
      ll_batched_device, ll_ptrs_device, info_batched, *magma_queue);

  zisa::copy(ll_batched_host, ll_batched_device);
  zisa::copy(ll_ptrs_host, ll_ptrs_device);

  for (int_t k = 0; k < n_mat; ++k) {
    for (int_t j = 0; j < n_cols; ++j) {
      for (int_t i = j; i < n_cols; ++i) {
        REQUIRE(
            zisa::almost_equal(llts[k](i, j), ll_batched_host(i, j, k), 1e-10));
      }
    }
  }
}

void check_random_matrices_against_eigen(int_t n_mat,
                                         int_t n_cols,
                                         int_t n_rhs) {
  auto matrices = make_spd_matrices(n_cols, n_mat);
  auto rhss = make_random_matrices(n_cols, n_rhs, n_mat);

  auto soln = std::vector<Eigen::MatrixXd>{};
  soln.reserve(n_mat);

  for (int_t k = 0; k < n_mat; ++k) {
    soln.emplace_back(matrices[k].llt().solve(rhss[k]));
  }

  auto ll_batched_host = convert_to_batched(matrices);
  auto ll_batched_device = array<double, 3, column_major>(
      ll_batched_host.shape(), device_type::cuda);
  zisa::copy(ll_batched_device, ll_batched_host);

  auto ll_ptrs_device = array<double *, 1>(n_mat, device_type::cuda);
  zisa::cuda::magma::fill_ptrs(ll_ptrs_device, ll_batched_device);
  auto info_batched = array<magma_int_t, 1>(n_mat, device_type::cuda);
  auto magma_queue = zisa::cuda::magma::make_default_queue();

  decompose_cholesky_inplace(
      ll_batched_device, ll_ptrs_device, info_batched, *magma_queue);

  auto rhs_batched_host = convert_to_batched(rhss);
  auto rhs_batched_device = array<double, 3, column_major>(
      rhs_batched_host.shape(), device_type::cuda);
  zisa::copy(rhs_batched_device, rhs_batched_host);

  auto rhs_ptrs_device = array<double *, 1>(n_mat, device_type::cuda);
  zisa::cuda::magma::fill_ptrs(rhs_ptrs_device, rhs_batched_device);

  solve_cholesky_inplace(ll_batched_device,
                         ll_ptrs_device,
                         rhs_batched_device,
                         rhs_ptrs_device,
                         *magma_queue);

  zisa::copy(rhs_batched_host, rhs_batched_device);

  for (int_t k = 0; k < n_mat; ++k) {
    for (int_t j = 0; j < n_rhs; ++j) {
      for (int_t i = 0; i < n_cols; ++i) {
        REQUIRE(
            zisa::almost_equal(soln[k](i, j), rhs_batched_host(i, j, k), 1e-8));
      }
    }
  }
}

}

TEST_CASE("CholeskyDecomposition; LLT vs Eigen", "[magma]") {
  zisa::int_t n_mat = 3;
  zisa::int_t n_cols = 5;

  zisa::check_llt_against_eigen(n_mat, n_cols);
}

TEST_CASE("CholeskyDecomposition; LLT solve vs Eigen", "[magma]") {

  zisa::int_t n_mat = 5;
  zisa::int_t n_cols = 3;
  zisa::int_t n_rhs = 7;

  zisa::check_random_matrices_against_eigen(n_mat, n_cols, n_rhs);
}
