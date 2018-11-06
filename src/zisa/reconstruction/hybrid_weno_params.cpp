#include <zisa/io/format_as_list.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>

namespace zisa {

HybridWENO_Params::HybridWENO_Params(StencilFamilyParams stencil_family_params,
                                     std::vector<double> linear_weights)
    : stencil_family_params(std::move(stencil_family_params)),
      linear_weights(std::move(linear_weights)) {}

std::ostream &operator<<(std::ostream &os, const HybridWENO_Params &params) {
  os << params.stencil_family_params << "\n";
  os << "linear weights: ";
  os << format_as_list(params.linear_weights);

  return os;
}

} // namespace zisa
