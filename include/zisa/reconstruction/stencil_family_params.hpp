#ifndef STENCIL_FAMILY_PARAMS_H_IN4IV
#define STENCIL_FAMILY_PARAMS_H_IN4IV

#include <cassert>
#include <string>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

struct StencilFamilyParams {
  std::vector<int> orders;
  std::vector<std::string> biases;
  std::vector<double> overfit_factors;

  int_t n_stencils() const;

public:
  StencilFamilyParams() = default;
  StencilFamilyParams(const StencilFamilyParams &) = default;
  StencilFamilyParams(StencilFamilyParams &&) = default;

  StencilFamilyParams(std::vector<int> orders,
                      std::vector<std::string> biases,
                      std::vector<double> overfit_factors);

  StencilFamilyParams &operator=(const StencilFamilyParams &) = default;
  StencilFamilyParams &operator=(StencilFamilyParams &&) = default;
};

StencilParams extract(const StencilFamilyParams &family_params, int_t k);

} // namespace zisa
#endif /* end of include guard */
