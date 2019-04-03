#ifndef ZISA_NEWTON_NKCOS_HPP
#define ZISA_NEWTON_NKCOS_HPP

#include <functional>

#include <zisa/math/basic_functions.hpp>

namespace zisa {

template <class F_DF>
double newton(const F_DF &f_df, double x_guess, double atol, int max_iter) {
  double x = x_guess;
  auto [fx, dfx] = f_df(x);

  int iter = 0;
  while (zisa::abs(fx / dfx) >= atol || iter >= max_iter) {
    // 0 ~ f(x) = f(x0) + df(x0)*(x - x0)
    // --> x = x0 - f(x0)/df(x0)

    x = x - fx / dfx;
    std::tie(fx, dfx) = f_df(x);

    ++iter;
  }

  LOG_ERR_IF(zisa::abs(fx / dfx) >= atol, "Failed to converge.");
  return x;
}

double newton(const std::function<double(double)> &f,
              const std::function<double(double)> &df,
              double x_guess,
              double atol = 1e-8,
              int max_iter = 100);

}

#endif // ZISA_NEWTON_HPP
