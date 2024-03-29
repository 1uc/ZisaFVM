// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <array>
#include <map>
#include <tuple>
#include <vector>

#include <zisa/math/triangular_rule.hpp>

namespace zisa {
static std::tuple<std::vector<double>, std::vector<std::array<double, 3>>>
permutate(double w, const std::vector<double> &lambda);

TriangularRule::TriangularRule(int_t n_points)
    : weights(shape_t<1>{n_points}), points(shape_t<1>{n_points}) {}

TriangularRule
make_quadrature_rule(const std::vector<double> &w,
                     const std::vector<std::array<double, 3>> &x) {
  assert(w.size() == x.size());

  auto qr = TriangularRule(w.size());

  for (int_t i = 0; i < w.size(); ++i) {
    qr.weights[i] = w[i];
    qr.points[i] = Barycentric2D{x[i][0], x[i][1], x[i][2]};
  }

  return qr;
}

TriangularRule make_triangular_rule(int_t deg) {
  if (deg == MAX_QUADRATURE_DEGREE) {
    deg = MAX_TRIANGULAR_RULE_DEGREE;
  }

  if (deg == 1) {
    auto [w, x] = permutate(1.0, {1.0 / 3.0});
    return make_quadrature_rule(w, x);

  } else if (deg == 2) {
    auto [w, x] = permutate(1.0 / 3.0, {2.0 / 3.0, 1.0 / 6.0});
    return make_quadrature_rule(w, x);

  } else if (deg == 3) {
    auto [w1, x1] = permutate(-0.5625, {1.0 / 3.0});
    auto [w2, x2] = permutate(1.5625 / 3.0, {0.6, 0.2});

    w1.reserve(w1.size() + w2.size());
    w1.insert(w1.end(), w2.begin(), w2.end());

    x1.reserve(x1.size() + x2.size());
    x1.insert(x1.end(), x2.begin(), x2.end());

    return make_quadrature_rule(w1, x1);
  } else if (deg == 4) {

    auto [w1, x1]
        = permutate(0.109951743655322, {0.816847572980459, 0.091576213509771});
    auto [w2, x2]
        = permutate(0.223381589678011, {0.108103018168070, 0.445948490915965});

    w1.reserve(w1.size() + w2.size());
    w1.insert(w1.end(), w2.begin(), w2.end());

    x1.reserve(x1.size() + x2.size());
    x1.insert(x1.end(), x2.begin(), x2.end());

    return make_quadrature_rule(w1, x1);
  } else if (deg == 5) {

    auto [w1, x1] = permutate(0.225, {1.0 / 3.0});
    auto [w2, x2]
        = permutate(0.132394152788506, {0.059715871789770, 0.470142064105115});
    auto [w3, x3]
        = permutate(0.125939180544827, {0.797426985353087, 0.101286507323456});

    w1.reserve(w1.size() + w2.size() + w3.size());
    w1.insert(w1.end(), w2.begin(), w2.end());
    w1.insert(w1.end(), w3.begin(), w3.end());

    x1.reserve(x1.size() + x2.size() + x3.size());
    x1.insert(x1.end(), x2.begin(), x2.end());
    x1.insert(x1.end(), x3.begin(), x3.end());

    return make_quadrature_rule(w1, x1);
  }

  LOG_ERR("Implement the missing case.");
}

const TriangularRule &cached_triangular_quadrature_rule(int_t deg) {
  static std::map<int_t, TriangularRule> qr;

  if (qr.find(deg) == qr.end()) {
    qr.insert({deg, make_triangular_rule(deg)});
  }

  return qr.at(deg);
}

static std::tuple<std::vector<double>, std::vector<std::array<double, 3>>>
permutate(double w, const std::vector<double> &lambda) {

  auto n = lambda.size();

  if (n == 1) {
    return {{w}, {{lambda[0], lambda[0], lambda[0]}}};
  } else if (n == 2) {
    return {{w, w, w},
            {{lambda[0], lambda[1], lambda[1]},
             {lambda[1], lambda[0], lambda[1]},
             {lambda[1], lambda[1], lambda[0]}}};
  } else if (n == 3) {
    return {{w, w, w, w, w, w},
            {{lambda[0], lambda[1], lambda[2]},
             {lambda[1], lambda[0], lambda[2]},

             {lambda[0], lambda[2], lambda[1]},
             {lambda[0], lambda[2], lambda[1]},

             {lambda[2], lambda[1], lambda[0]},
             {lambda[1], lambda[2], lambda[0]}}};
  }

  LOG_ERR("Implement the missing case.");
}

} // namespace zisa
