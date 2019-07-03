#ifndef WENO_AO_PARAMS_H_SH3SV
#define WENO_AO_PARAMS_H_SH3SV

#include <vector>

#include <zisa/config.hpp>
#include <zisa/reconstruction/stencil_family_params.hpp>

namespace zisa {

struct HybridWENOParams {

  HybridWENOParams() = default;
  HybridWENOParams(const HybridWENOParams &) = default;
  HybridWENOParams(HybridWENOParams &&) = default;

  HybridWENOParams(StencilFamilyParams stencil_family_params,
                   std::vector<double> linear_weights,
                   double epsilon,
                   double exponent);

  HybridWENOParams &operator=(const HybridWENOParams &) = default;
  HybridWENOParams &operator=(HybridWENOParams &&) = default;

public:
  StencilFamilyParams stencil_family_params;
  std::vector<double> linear_weights;
  double epsilon;
  double exponent;
};

int max_order(const HybridWENOParams &params);

HybridWENOParams make_hybrid_weno_params(int order);

std::ostream &operator<<(std::ostream &os, const HybridWENOParams &params);

} // namespace zisa
#endif /* end of include guard */
