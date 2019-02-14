#ifndef TRIANGULAR_RULE_H_V7Y5S
#define TRIANGULAR_RULE_H_V7Y5S

#include <zisa/config.hpp>
#include <zisa/math/barycentric.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

struct TriangularRule {
  array<double, 1> weights;
  array<Barycentric, 1> points;

  TriangularRule(int_t n_points);
  TriangularRule(const TriangularRule &qr) = default;
  TriangularRule(TriangularRule &&qr) = default;
};

TriangularRule make_triangular_rule(int_t deg);

const TriangularRule &cached_triangular_quadrature_rule(int_t deg);

} // namespace zisa
#endif
