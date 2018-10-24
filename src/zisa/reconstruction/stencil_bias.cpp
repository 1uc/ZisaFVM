#include <zisa/utils/string_format.hpp>
#include <zisa/utils/logging.hpp>
#include <zisa/reconstruction/stencil_bias.hpp>

namespace zisa {

StencilBias deduce_bias(const std::string &b) {

  if (b.compare("c") == 0) {
    return StencilBias::central;
  } else if (b.compare("b") == 0) {
    return StencilBias::one_sided;
  }

  LOG_ERR(string_format("Missing case. [%s]", b.c_str()));
}

} // namespace zisa
