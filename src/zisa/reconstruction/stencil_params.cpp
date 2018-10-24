#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

StencilParams::StencilParams(int order, std::string bias, double overfit_factor)
    : order(order), bias(bias), overfit_factor(overfit_factor) {}

} // namespace zisa
