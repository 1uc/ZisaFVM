// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef STENCIL_H_FGW9W
#define STENCIL_H_FGW9W

#include <vector>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/cone.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/reconstruction/stencil_bias.hpp>
#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

class Stencil {
public:
  Stencil() = default;

  /// Simplified constructor for one point stencils.
  explicit Stencil(int_t i_cell);

  Stencil(std::vector<int_t> &global_indices, const StencilParams &params);

  Stencil(std::vector<int_t> &l2g,
          const Grid &grid,
          int_t i_cell,
          const StencilParams &params);

  Stencil(std::vector<int_t> &l2g,
          const Grid &grid,
          int_t i_cell,
          int_t k,
          const StencilParams &params);

  int_t local(int_t k) const;
  int_t global(int_t k) const;

  const array<int_t, 1> &local() const { return local_; }
  const array<int_t, 1> &global() const { return global_; }

  /// Factor by which the LSQ problem is over determined.
  double overfit_factor() const;

  /// What is the bias of this stencil (central, one-sided)?
  StencilBias bias() const;

  /// Attainable order with this stencil.
  int order() const;

  /// Requested (maximum) order for this stencil.
  int max_order() const;

  /// Size of the stencil.
  int_t size() const;

  /// Apply a permutation to the global indices.
  template <class F>
  void apply_permutation(const F &f) {
    for (auto &i : global_) {
      i = f(i);
    }
  }

protected:
  int_t max_size() const;
  void assign_local_indices(const std::vector<int_t> &global_indices,
                            std::vector<int_t> &l2g);

private:
  int max_order_;
  int order_;

  int_t size_;
  int_t max_size_;

  StencilBias bias_;
  double overfit_factor_;

  array<int_t, 1> local_;
  array<int_t, 1> global_;
};

bool operator==(const Stencil &a, const Stencil &b);
bool operator!=(const Stencil &a, const Stencil &b);

std::vector<int_t> central_stencil(const Grid &grid, int_t i, int_t n_points);

/// A biased stencil which will result in a valid LSQ matrix.
std::vector<int_t> biased_stencil(
    const Grid &grid, int_t i_center, int_t k, int_t n_points, int order);

/// Biased stencil using condidates in a wider cone.
std::vector<int_t> less_conservative_stencil(const Grid &grid,
                                             int_t i_center,
                                             int_t k,
                                             int_t n_points);

/// Biased stencil using only candidates form the ideal cone.
std::vector<int_t>
conservative_stencil(const Grid &grid, int_t i_center, int_t k, int_t n_points);

/// A stencil consisting of points intersecting the `region`.
std::vector<int_t> region_based_stencil(const Grid &grid,
                                        int_t i,
                                        int_t n_points,
                                        const Region &region);

int deduce_max_order(int_t stencil_size, double factor, int n_dims);
int_t required_stencil_size(int deg, double factor, int n_dims);

} // namespace zisa

#endif /* end of include guard */
