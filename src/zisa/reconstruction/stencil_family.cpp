#include <limits>

#include <zisa/math/comparison.hpp>
#include <zisa/reconstruction/stencil.hpp>
#include <zisa/reconstruction/stencil_family.hpp>

namespace zisa {

StencilFamily::StencilFamily(const std::shared_ptr<Grid> &grid,
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

  combined_stencil_size_ = l2g.size();

  order_ = 0;
  for (int_t i = 0; i < size(); ++i) {
    order_ = zisa::max(order_, (*this)[i].order());
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

bool operator==(const StencilFamily &lhs, const StencilFamily &rhs) {
  if (lhs.combined_stencil_size() != rhs.combined_stencil_size()) {
    return false;
  }

  if (!std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end())) {
    return false;
  }

  auto l2g_lhs = lhs.local2global();
  auto l2g_rhs = rhs.local2global();

  auto stencil_size = lhs.combined_stencil_size();

  return std::equal(lhs.begin(),
                    lhs.begin() + stencil_size,
                    rhs.begin(),
                    rhs.begin() + stencil_size);
}

bool operator!=(const StencilFamily &lhs, const StencilFamily &rhs) {
  return !(lhs == rhs);
}

int_t highest_order_central_stencil(const StencilFamily &stencils) {
  int_t k_high = 0;

  for (int_t i = 1; i < stencils.size(); ++i) {
    if (stencils[i].order() > stencils[k_high].order()
        && stencils[i].bias() == StencilBias::central) {
      k_high = i;
    }
  }

  assert(stencils[k_high].bias() == StencilBias::central);
  return k_high;
}

} // namespace zisa