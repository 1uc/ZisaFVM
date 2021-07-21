// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_TETRAHEDRAL_RULE_HPP_CHUIWKJ
#define ZISA_TETRAHEDRAL_RULE_HPP_CHUIWKJ

#include <zisa/config.hpp>

#include <zisa/math/barycentric.hpp>
#include <zisa/math/max_quadrature_degree.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {
struct TetrahedralRule {
  array<double, 1> weights;
  array<Barycentric3D, 1> points;

  TetrahedralRule() = default;
  TetrahedralRule(const TetrahedralRule &other) = default;
  TetrahedralRule(TetrahedralRule &&other) noexcept = default;

  explicit TetrahedralRule(int_t n_points)
      : weights(shape_t<1>{n_points}), points(shape_t<1>{n_points}) {}

  TetrahedralRule(array<double, 1> weights, array<Barycentric3D, 1> &points)
      : weights(std::move(weights)), points(std::move(points)) {}

  TetrahedralRule &operator=(const TetrahedralRule &other) = default;
  TetrahedralRule &operator=(TetrahedralRule &&other) noexcept = default;
};

/// Construct a quadrature rule on a tetrahedron of order at least `deg`.
/** This generates optimal, symmetric quadrature rules with positive weights.
 *
 *  Reference:
 *    L. Shunn, F. Ham, Journal of Computational and Applied Mathematics, 2012.
 */
TetrahedralRule make_tetrahedral_rule(int_t deg);

const TetrahedralRule &cached_tetrahedral_rule(int_t deg);

constexpr int_t MAX_TETRAHEDRAL_RULE_DEGREE = 3;

}

#endif // ZISA_TETRAHEDRAL_RULE_HPP
