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
};

const EdgeRule &cached_edge_quadrature_rule(int_t deg);

} // namespace zisa

#endif
