// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include "zisa/math/triangular_rule.hpp"
#include <zisa/config.hpp>
#include <zisa/cuda/hello_world.hpp>
#include <zisa/cuda/reconstruction/qr_decomposition.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/max_quadrature_degree.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/assemble_weno_ao_matrix.hpp>
#include <zisa/reconstruction/lsq_solver.hpp>
#include <zisa/reconstruction/qr_decomposition.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

#include <Eigen/Dense>
#include <algorithm>
#include <magma_v2.h>
#include <map>
#include <vector>

namespace zisa {

bool contains(const std::vector<int_t> &cells, int_t i) {
  return std::find(cells.begin(), cells.end(), i) != cells.end();
}

bool contains(const std::map<std::pair<int_t, int_t>, int_t> &m,
              const std::pair<int_t, int_t> &k) {
  return m.find(k) != m.end();
}

void test_batched_qr() {
  auto grid = zisa::load_grid(
      "/scratch/lucg/zisa-grids/grids/circle-01/00/grid.msh.h5",
      MAX_TRIANGULAR_RULE_DEGREE);

  auto n_dims = grid->n_dims();

  auto params
      = StencilFamilyParams{{3, 2, 2}, {"c", "b", "b"}, {2.0, 1.5, 1.5}};

  auto stencil_families = zisa::compute_stencil_families(*grid, params);

  std::map<std::pair<int_t, int_t>, int_t> size_to_bin;
  std::vector<std::vector<std::pair<int_t, int_t>>> stencil_indices;
  std::vector<std::vector<Stencil>> binned_stencils;

  for (auto i : cell_indices(*grid)) {
    for (auto k : index_range(stencil_families[i].size())) {
      const auto &stencil = stencil_families[i][k];

      auto bin_key = std::pair<int_t, int_t>{stencil.size(), stencil.order()};

      if (!contains(size_to_bin, bin_key)) {
        size_to_bin[bin_key] = size_to_bin.size();
        stencil_indices.emplace_back();
        binned_stencils.emplace_back();
      }

      auto bin = size_to_bin[bin_key];

      binned_stencils[bin].push_back(stencil);
      stencil_indices[bin].emplace_back(i, k);
    }
  }

  for (auto k : size_to_bin) {
    PRINT(k.first.first);
    PRINT(k.first.second);
    PRINT(k.second);
    PRINT(stencil_indices[k.second].size());
    PRINT(binned_stencils[k.second].size());
  }

  std::vector<array<double, 3>> binned_weno_matrices(size_to_bin.size());
  for (auto kv : size_to_bin) {
    auto bin = kv.second;
    auto n_mat = binned_stencils[bin].size();

    auto n_points = binned_stencils[bin][0].size();
    auto order = binned_stencils[bin][0].order();

    auto degree = order - 1;
    auto n_rows = n_points - 1;
    auto n_cols = poly_dof(degree, n_dims) - 1;

    if (order == 1) {
      binned_weno_matrices[bin] = array<double, 3>(shape_t<3>{1, 1, 1});
      binned_weno_matrices[bin](0, 0, 0) = 1.0;
    } else {
      binned_weno_matrices[bin]
          = array<double, 3>(shape_t<3>{n_mat, n_rows, n_cols});

      for (auto k_mat : index_range(binned_stencils.size())) {
        const auto stencil = binned_stencils[bin][k_mat];

        LOG_ERR_IF(stencil.size() != n_points,
                   "Failed to bin stencils properly by size.");

        LOG_ERR_IF(stencil.order() != order,
                   "Failed to bin stencils properly by order.");

        auto raw_ptr = &binned_weno_matrices[bin](k_mat, 0, 0);
        Eigen::Map<Eigen::MatrixXd> A(raw_ptr, n_rows, n_cols);
        Eigen::Ref<Eigen::MatrixXd> A_ref(A);
        assemble_weno_ao_matrix(A_ref, *grid, stencil.global(), order);
      }
    }
  }

  auto d_binned_weno_matrices = std::vector<array<double, 3>>();
  for (auto k : index_range(binned_weno_matrices.size())) {
    d_binned_weno_matrices.emplace_back(binned_weno_matrices[k].shape(),
                                        device_type::cuda);

    zisa::copy(d_binned_weno_matrices[k], binned_weno_matrices[k]);
  }

  auto magma_context
      = std::make_shared<zisa::cuda::MAGMAQueue>(magma_int_t(0));

  auto qr_decomposition
      = zisa::cuda::QRDecomposition(magma_context, d_binned_weno_matrices[0]);
}
}

int main() {
#if ZISA_HAS_CUDA > 0
  magma_init();

  zisa::hello_world();
  zisa::test_batched_qr();

  magma_finalize();

#else
  LOG_ERR("Compiled without CUDA support.");
#endif
}
