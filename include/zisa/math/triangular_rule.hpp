#ifndef TRIANGULAR_RULE_H_V7Y5S
#define TRIANGULAR_RULE_H_V7Y5S

#include <zisa/config.hpp>
#include <zisa/math/barycentric.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

struct TriangularRule {
  array<double, 1> weights;
  array<Barycentric2D, 1> points;

  TriangularRule(int_t n_points);
  TriangularRule(const TriangularRule &qr) = default;
  TriangularRule(TriangularRule &&qr) = default;
};

/// Compute weights and quadrature points.
/**
 * References:
 *   [1] D.A. Dunavant, High degree efficient symmetrical Gaussian quadrature
 *       rules for the triangle, 1985.
 */
TriangularRule make_triangular_rule(int_t deg);

constexpr int_t MAX_TRIANGULAR_RULE_DEGREE = 5;

const TriangularRule &cached_triangular_quadrature_rule(int_t deg);

} // namespace zisa
#endif
