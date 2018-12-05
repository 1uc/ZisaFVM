#include <zisa/math/edge_rule.hpp>
#include <zisa/math/gauss_legendre.hpp>

namespace zisa {

EdgeRule::EdgeRule(int_t deg)

    : weights((deg + 1) / 2), points((deg + 1) / 2) {

  int_t n_points = (deg + 1) / 2;

  if (n_points == 1) {
    auto qr = make_gauss_legendre<1>();

    for (int_t k = 0; k < n_points; ++k) {
      weights[k] = 0.5 * qr.weights[k];
      points[k] = qr.points[k];
    }
  } else if (n_points == 2) {
    auto qr = make_gauss_legendre<2>();

    for (int_t k = 0; k < n_points; ++k) {
      weights[k] = 0.5 * qr.weights[k];
      points[k] = qr.points[k];
    }
  } else if (n_points == 3) {
    auto qr = make_gauss_legendre<3>();

    for (int_t k = 0; k < n_points; ++k) {
      weights[k] = 0.5 * qr.weights[k];
      points[k] = qr.points[k];
    }
  } else if (n_points == 4) {
    auto qr = make_gauss_legendre<4>();

    for (int_t k = 0; k < n_points; ++k) {
      weights[k] = 0.5 * qr.weights[k];
      points[k] = qr.points[k];
    }
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
