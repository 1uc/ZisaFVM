#include <limits>

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

  order_ = 0;
  for (int_t i = 0; i < size(); ++i) {
    auto stencil_order = (*this)[i].order();

//    if(stencil_order == 1 && params.orders[i] != 1) {
//      PRINT("truncating");
//      truncate_all_stencils_to_first_order(i_cell);
//      break;
//    }

    order_ = std::max(order_, stencil_order);
  }

  combined_stencil_size_ = l2g.size();
  k_high_ = highest_order_central_stencil(*this);
}

void StencilFamily::truncate_all_stencils_to_first_order(int_t i_cell) {
  for(auto &stencil : stencils_) {
    stencil = Stencil(i_cell);
  }
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

auto StencilFamily::begin() const -> decltype(stencils_.begin()) {
  return stencils_.begin();
}
auto StencilFamily::end() const -> decltype(stencils_.end()) {
  return stencils_.end();
}

int_t StencilFamily::highest_order_stencil() const { return k_high_; }

bool operator==(const StencilFamily &lhs, const StencilFamily &rhs) {
  if (lhs.combined_stencil_size() != rhs.combined_stencil_size()) {
    return false;
  }

  return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

bool operator!=(const StencilFamily &lhs, const StencilFamily &rhs) {
  return !(lhs == rhs);
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
