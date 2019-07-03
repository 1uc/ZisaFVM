#ifndef ZISA_TETRAHEDRAL_RULE_HPP_CHUIWKJ
#define ZISA_TETRAHEDRAL_RULE_HPP_CHUIWKJ

#include <zisa/config.hpp>

#include <zisa/math/barycentric.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {
struct TetrahedralRule {
  array<double, 1> weights;
  array<Barycentric3D, 1> points;

  TetrahedralRule(int_t n_points)
      : weights(shape_t<1>{n_points}), points(shape_t<1>{n_points}) {}
};

/// Construct a quadrature rule on a tetrahedron of order at least `deg`.
/** This generates optimal, symmetric quadrature rules with positive weights.
 *
 *  Reference:
 *    L. Shunn, F. Ham, Journal of Computational and Applied Mathematics, 2012.
 */
TetrahedralRule make_tetrahedral_rule(int_t deg);

}

#endif // ZISA_TETRAHEDRAL_RULE_HPP
