#ifndef QUASI_NEWTON_H_W8GTO
#define QUASI_NEWTON_H_W8GTO

#include <zisa/math/cartesian_expr.hpp>
#include <zisa/math/rolling_convergence_rate.hpp>

namespace zisa {

template <class F, class IDF, class X>
std::tuple<X, bool> quasi_newton(const F &f,
                                 const IDF inv_df,
                                 const X &x0,
                                 const X &atol,
                                 int max_iter = 20) {

  X x = x0;
  auto fx = f(x0);

  X dx = X(atol + 1.0);

  auto is_converged = [&dx, &atol](double factor = 1.0) {
    return zisa::all(zisa::abs(dx) <= factor * atol);
  };

  auto rate = RollingConvergenceRate<X>(atol);

  int iter = 0;
  while (!is_converged() && iter < max_iter) {
    auto inv_df_x = inv_df(x);
    dx = inv_df_x(fx);
    x -= dx;

    rate.push_back(dx);

    if (iter >= 4 && !rate.is_converging(X(0.0))) {
      LOG_WARN("Is not converging.");
      PRINT(x);
      PRINT(dx);
      PRINT(fx);

      return std::tuple<X, bool>{x0, false};
    }

    fx = f(x);
    ++iter;
  }

  if (iter == max_iter && !is_converged(1000.0)) {
    LOG_WARN("Didn't converge.");
    return std::tuple<X, bool>{x0, false};
  }

  return std::tuple<X, bool>{x, true};
}

} // namespace zisa
#endif /* end of include guard */
