#ifndef WENO_AO_PARAMS_H_SH3SV
#define WENO_AO_PARAMS_H_SH3SV

#include <vector>

#include <zisa/config.hpp>
#include <zisa/reconstruction/stencil_family_params.hpp>

namespace zisa {

struct HybridWENO_Params {

  HybridWENO_Params() = default;
  HybridWENO_Params(const HybridWENO_Params &) = default;
  HybridWENO_Params(HybridWENO_Params &&) = default;

  HybridWENO_Params(StencilFamilyParams stencil_family_params,
                 std::vector<double> linear_weights);

  HybridWENO_Params &operator=(const HybridWENO_Params &) = default;
  HybridWENO_Params &operator=(HybridWENO_Params &&) = default;

public:
  StencilFamilyParams stencil_family_params;
  std::vector<double> linear_weights;
};

std::ostream& operator<<(std::ostream &os, const HybridWENO_Params &params);

} // namespace zisa
#endif /* end of include guard */
