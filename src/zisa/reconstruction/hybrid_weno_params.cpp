#include <zisa/io/format_as_list.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>

namespace zisa {

HybridWENOParams::HybridWENOParams(StencilFamilyParams stencil_family_params,
                                   std::vector<double> linear_weights,
                                   double epsilon,
                                   double exponent)
    : stencil_family_params(std::move(stencil_family_params)),
      linear_weights(std::move(linear_weights)),
      epsilon(epsilon),
      exponent(exponent) {}

std::ostream &operator<<(std::ostream &os, const HybridWENOParams &params) {
  os << params.stencil_family_params << "\n";
  os << "linear weights: ";
  os << format_as_list(params.linear_weights);

  return os;
}

HybridWENOParams make_hybrid_weno_params(int order) {
  if (order == 1) {
    return HybridWENOParams({{1}, {"c"}, {2.0}}, {1.0}, 1e-6, 4);
  }

  if (order == 2) {
    return HybridWENOParams(
        {{order, 2, 2, 2}, {"c", "b", "b", "b"}, {3.0, 2.0, 2.0, 2.0}},
        {100.0, 1.0, 1.0, 1.0},
        1e-6,
        4);
  }

  return HybridWENOParams(
      {{order, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
      {100.0, 1.0, 1.0, 1.0},
      1e-6,
      4);
}

} // namespace zisa
