#include <zisa/io/format_as_list.hpp>
#include <zisa/reconstruction/stencil_family_params.hpp>
#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

StencilFamilyParams::StencilFamilyParams(std::vector<int> orders,
                                         std::vector<std::string> biases,
                                         std::vector<double> overfit_factors)
    : orders(std::move(orders)),
      biases(std::move(biases)),
      overfit_factors(std::move(overfit_factors)) {

  auto size_ = this->orders.size();
  assert(size_ == this->biases.size());
  assert(size_ == this->overfit_factors.size());

  ZISA_UNUSED(size_);
}

int_t StencilFamilyParams::n_stencils() const {
  return orders.size();
}

StencilParams extract(const StencilFamilyParams &family_params, int_t k) {
  return {family_params.orders[k],
          family_params.biases[k],
          family_params.overfit_factors[k]};
}

std::ostream &operator<<(std::ostream &os, const StencilFamilyParams &params) {
  os << "orders: " << format_as_list(params.orders) << "\n";
  os << "biases: " << format_as_list(params.biases) << "\n";
  os << "overfit_factors: " << format_as_list(params.overfit_factors);

  return os;
}

int max_order(const StencilFamilyParams &params) {
  return *std::max_element(params.orders.begin(), params.orders.end());
}

} // namespace zisa
