/* Quadrature rules.
 */

#include <vector>
#include <zisa/math/quadrature.hpp>

namespace zisa {

std::tuple<std::vector<double>, std::vector<std::vector<double>>>
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

Barycentric::Barycentric(double lambda1, double lambda2, double lambda3)
    : lambda{lambda1, lambda2, lambda3} {}

std::ostream &operator<<(std::ostream &os, const Barycentric &bc) {
  os << string_format("[ %e, %e, %e]", bc[0], bc[1], bc[2]);
  return os;
}

XY Barycentric::operator()(const Triangle &tri) const {
  return XY(tri.A * lambda[0] + tri.B * lambda[1] + tri.C * lambda[2]);
}

QuadratureRule::QuadratureRule(int_t n_points)
    : weights(shape_t<1>{n_points}), points(shape_t<1>{n_points}) {}

QuadratureRule make_quadrature_rule(const std::vector<double> &w,
                                    const std::vector<std::vector<double>> &x) {
  assert(w.size() == x.size());

  auto qr = QuadratureRule(w.size());

  for (int_t i = 0; i < w.size(); ++i) {
    qr.weights[i] = w[i];
    qr.points[i] = Barycentric{x[i][0], x[i][1], x[i][2]};
  }

  return qr;
}

QuadratureRule make_triangular_rule(int deg) {
  if (deg == 1) {
    auto [w, x] = permutate(1.0, {1.0 / 3.0});
    return make_quadrature_rule(std::move(w), std::move(x));

  } else if (deg == 2) {
    auto [w, x] = permutate(1.0 / 3.0, {2.0 / 3.0, 1.0 / 6.0});
    return make_quadrature_rule(std::move(w), std::move(x));

  } else if (deg == 3) {
    auto [w1, x1] = permutate(-0.5625, {1.0 / 3.0});
    auto [w2, x2] = permutate(1.5625 / 3.0, {0.6, 0.2});

    w1.reserve(w1.size() + w2.size());
    w1.insert(w1.end(), w2.begin(), w2.end());

    x1.reserve(x1.size() + x2.size());
    x1.insert(x1.end(), x2.begin(), x2.end());

    return make_quadrature_rule(std::move(w1), std::move(x1));
  } else if (deg == 4) {

    auto [w1, x1]
        = permutate(0.109951743655322, {0.816847572980459, 0.091576213509771});
    auto [w2, x2]
        = permutate(0.223381589678011, {0.108103018168070, 0.445948490915965});

    w1.reserve(w1.size() + w2.size());
    w1.insert(w1.end(), w2.begin(), w2.end());

    x1.reserve(x1.size() + x2.size());
    x1.insert(x1.end(), x2.begin(), x2.end());

    return make_quadrature_rule(std::move(w1), std::move(x1));
  }

  LOG_ERR("Implement the missing case.");
}

} // namespace zisa
