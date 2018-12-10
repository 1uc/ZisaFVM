#ifndef EDGE_RULE_H_EFB8F
#define EDGE_RULE_H_EFB8F

#include<map>

#include <zisa/config.hpp>
#include <zisa/math/gauss_legendre.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

class EdgeRule {
public:
  array<double, 1> weights;
  array<double, 1> points;

public:
  EdgeRule(int_t deg);

  std::string str() const;

private:
  template <int_t deg>
  void init();
};

template <int_t n_points>
void EdgeRule::init() {
  auto qr = make_gauss_legendre<n_points>();

  for (int_t k = 0; k < n_points; ++k) {
    weights[k] = 0.5 * qr.weights[k];
    points[k] = qr.points[k];
  }
}

const EdgeRule &cached_edge_quadrature_rule(int_t deg);

} // namespace zisa

#endif
