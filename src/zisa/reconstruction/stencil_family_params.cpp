#include <zisa/reconstruction/stencil_family_params.hpp>
#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

StencilFamilyParams::StencilFamilyParams(std::vector<int> orders,
                                         std::vector<std::string> biases,
                                         std::vector<double> overfit_factors)
    : orders(std::move(orders)),
      biases(std::move(biases)),
      overfit_factors(std::move(overfit_factors)) {}

int_t StencilFamilyParams::n_stencils() const {
  assert(orders.size() == biases.size());
  assert(overfit_factors.size() == biases.size());

  return orders.size();
}

StencilParams extract(const StencilFamilyParams &family_params, int_t k) {
  return {family_params.orders[k],
          family_params.biases[k],
          family_params.overfit_factors[k]};
}

} // namespace zisa
