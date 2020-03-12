#include <zisa/math/edge_rule.hpp>
#include <zisa/math/gauss_legendre.hpp>

namespace zisa {

EdgeRule::EdgeRule(int_t deg)

    : weights(deg / 2 + 1), points(deg / 2 + 1) {

  int_t n_points = deg / 2 + 1;

  if (n_points == 1) {
    init<1>();
  } else if (n_points == 2) {
    init<2>();
  } else if (n_points == 3) {
    init<3>();
  } else if (n_points == 4) {
    init<4>();
  } else {
    LOG_ERR(string_format("Implement case. [%d]", deg));
  }
}

std::string EdgeRule::str() const {
  auto n_points = weights.size();
  return string_format("Edge quadrature: %d point rule", n_points);
}

const EdgeRule &cached_edge_quadrature_rule(int_t deg) {
  static std::map<int_t, EdgeRule> qr;

  if (qr.find(deg) == qr.end()) {
    qr.insert({deg, EdgeRule(deg)});
  }

  return qr.at(deg);
}

} // namespace zisa
