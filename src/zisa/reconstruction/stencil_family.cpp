// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <limits>

#include <zisa/loops/for_each.hpp>
#include <zisa/math/comparison.hpp>
#include <zisa/reconstruction/stencil.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {

int_t highest_order_central_stencil(const StencilFamily &stencils);

StencilFamily::StencilFamily(const Grid &grid,
                             int_t i_cell,
                             const StencilFamilyParams &params)
    : stencils_(shape_t<1>{params.n_stencils()}) {

  int_t k_biased = 0;
  for (int_t i = 0; i < size(); ++i) {
    if (deduce_bias(params.biases[i]) == StencilBias::one_sided) {
      stencils_(i) = Stencil(l2g, grid, i_cell, k_biased, extract(params, i));
      ++k_biased;
    } else {
      stencils_(i) = Stencil(l2g, grid, i_cell, extract(params, i));
    }
  }

  //  auto d = distance_to_boundary(grid, i_cell, /* max_distance = */ 2);
  //  if (d <= 1) {
  //    truncate_all_stencils_to_first_order(i_cell);
  //  }

  order_ = 1;
  for (int_t i = 0; i < size(); ++i) {
    auto stencil_order = (*this)[i].order();
    order_ = std::max(order_, stencil_order);
  }

  combined_stencil_size_ = l2g.size();
  k_high_ = highest_order_central_stencil(*this);
}

const Stencil &StencilFamily::operator[](int_t k) const {
  assert(k < size());
  return stencils_[k];
}

const std::vector<int_t> &StencilFamily::local2global() const { return l2g; }

int_t StencilFamily::size() const { return stencils_.size(); }

int_t StencilFamily::combined_stencil_size() const {
  return combined_stencil_size_;
}

int StencilFamily::order() const { return order_; }

auto StencilFamily::begin() -> decltype(stencils_.begin()) {
  return stencils_.begin();
}

auto StencilFamily::begin() const -> decltype(stencils_.begin()) {
  return stencils_.begin();
}

auto StencilFamily::end() -> decltype(stencils_.end()) {
  return stencils_.end();
}

auto StencilFamily::end() const -> decltype(stencils_.end()) {
  return stencils_.end();
}

int_t StencilFamily::highest_order_stencil() const { return k_high_; }

void StencilFamily::truncate_to_first_order(int_t i_cell) {
  stencils_ = array<Stencil, 1>(1);
  stencils_[0] = Stencil(i_cell);
  l2g = std::vector<int_t>{i_cell};
  combined_stencil_size_ = 1;
  order_ = 1;
  k_high_ = 0;
}

bool operator==(const StencilFamily &lhs, const StencilFamily &rhs) {
  if (lhs.combined_stencil_size() != rhs.combined_stencil_size()) {
    return false;
  }

  return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

bool operator!=(const StencilFamily &lhs, const StencilFamily &rhs) {
  return !(lhs == rhs);
}

array<StencilFamily, 1>
compute_stencil_families(const Grid &grid, const StencilFamilyParams &params) {
  auto stencil_families = array<StencilFamily, 1>(grid.n_cells);
  // Empirically, this loop *must* be serial or else the code will segfault,
  // occasionally.
  // TODO Find out why, then make it parallel.
  for_each(serial_policy{},
           flat_range(stencil_families),
           [&stencil_families, &grid, &params](int_t i) {
             if (grid.cell_flags[i].interior
                 || grid.cell_flags[i].ghost_cell_l1) {
               stencil_families(i) = StencilFamily(grid, i, params);
             } else {
               auto o1 = StencilFamilyParams{{1}, {"c"}, {1.0}};
               stencil_families(i) = StencilFamily(grid, i, o1);
             }
           });

  return stencil_families;
}

int_t highest_order_central_stencil(const StencilFamily &stencils) {
  int_t k_high = 0;

  for (int_t i = 1; i < stencils.size(); ++i) {

    auto k_order = stencils[k_high].order();

    if (stencils[i].order() > k_order) {
      k_high = i;
    } else if (stencils[i].order() == k_order
               && stencils[i].bias() == StencilBias::central) {
      k_high = i;
    }
  }

  return k_high;
}

} // namespace zisa
