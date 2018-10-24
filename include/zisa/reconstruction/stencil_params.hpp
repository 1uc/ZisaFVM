#ifndef STENCIL_PARAMS_H_UXTH0
#define STENCIL_PARAMS_H_UXTH0

#include <string>

namespace zisa {

struct StencilParams {
  StencilParams() = default;
  StencilParams(const StencilParams &) = default;
  StencilParams(StencilParams &&) = default;

  StencilParams(int order, std::string bias, double overfit_factor);

  StencilParams &operator=(const StencilParams &) = default;
  StencilParams &operator=(StencilParams &&) = default;

  int order;
  std::string bias;
  double overfit_factor;
};

} // namespace zisa
#endif /* end of include guard */
