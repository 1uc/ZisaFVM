/*
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2017-07-24
 */

#ifndef FOURIER_SERIES_H_R3XD1
#define FOURIER_SERIES_H_R3XD1

#include <zisa/config.hpp>
#include <zisa/math/polynomial_expr.hpp>

namespace zisa {

template <int DEGREE> class FourierSeriesBase {
public:
  static constexpr int degree() { return DEGREE; }
  static constexpr int size() { return degree() + 1; }

public:
  FourierSeriesBase(void) = default;

  FourierSeriesBase(const std::initializer_list<double> &init_list) {
    assert(int(init_list.size()) == size());
    std::copy(init_list.begin(), init_list.end(), coeffs);
  }

  template <class E> FourierSeriesBase(const PolynomialCRTP<E> &e_) {
    deep_copy(static_cast<const E &>(e_));
  }

  template <class E> void operator=(const PolynomialCRTP<E> &e_) {
    const E &e = static_cast<const E &>(e_);

    deep_copy(e);
  }

  template <class E> void deep_copy(const PolynomialCRTP<E> &e_) {
    const E &e = static_cast<const E &>(e_);
    static_assert(E::degree() <= degree(), "Polynomial degree mismatch.");

    for (int i = 0; i <= degree(); ++i) {
      this->at(i) = e.at(i);
    }
  }

  double &at(int i) { return a(i); }

  double &a(int i) { return coeffs[i]; }

  double at(int i) const { return a(i); }

  double a(int i) const {
    assert(i >= 0);

    if (i <= degree()) {
      return coeffs[i];
    } else {
      return 0.0;
    }
  }

protected:
  double coeffs[size()];
};

template <int DEGREE> class FourierSeriesDerivative;

template <int DEGREE>
class FourierSeries : public FourierSeriesBase<DEGREE>,
                      public PolynomialCRTP<FourierSeries<DEGREE>> {
private:
  using super = FourierSeriesBase<DEGREE>;

public:
  using super::FourierSeriesBase;

  double operator()(double alpha) const {
    double fx = this->a(0); // cos(0) = 1

    for (int i = 1; i <= this->degree(); ++i) {
      fx += this->a(i) * zisa::cos(i * alpha);
    }

    return fx;
  }

  FourierSeriesDerivative<DEGREE> derivative(void) const {
    return FourierSeriesDerivative<DEGREE>(*this);
  }
};

template <int DEGREE>
class FourierSeriesDerivative : public FourierSeriesBase<DEGREE>,
                                public PolynomialCRTP<FourierSeries<DEGREE>> {
private:
  using super = FourierSeriesBase<DEGREE>;

public:
  using super::FourierSeriesBase;

  double operator()(double alpha) const {
    double fx = 0.0;

    for (int i = 1; i <= this->degree(); ++i) {
      fx += -this->a(i) * i * zisa::sin(i * alpha);
    }

    return fx;
  }

  FourierSeriesDerivative derivative(void) const {
    return FourierSeriesDerivative(*this);
  }
};

} // namespace zisa
#endif /* end of include guard */
