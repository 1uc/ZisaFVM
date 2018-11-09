/* Simple arrays.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2014-09-15
 */
#ifndef CARTESIAN_EXPR_H_KJMKGEYM
#define CARTESIAN_EXPR_H_KJMKGEYM

#include "zisa/config.hpp"
#include <cmath>

namespace zisa {

// CRTP Base
template <class E, class T>
class CartesianExpr {
public:
  using scalar_t = T;
};

template <class E1, class E2, class T>
class CartesianPlus : public CartesianExpr<CartesianPlus<E1, E2, T>, T> {
public:
  using scalar_t = T;

  ANY_DEVICE_INLINE CartesianPlus(const CartesianExpr<E1, T> &e1,
                                  const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  ANY_DEVICE_INLINE scalar_t operator()(int_t i) const { return e1(i) + e2(i); }

  ANY_DEVICE_INLINE constexpr static int_t size(void) {
    static_assert(E1::size() == E2::size(), "Length mismatch.");
    return E1::size();
  }

private:
  const E1 &e1;
  const E2 &e2;
};

template <class E1, class E2, class T>
class CartesianProduct : public CartesianExpr<CartesianProduct<E1, E2, T>, T> {
public:
  using scalar_t = T;

  ANY_DEVICE_INLINE CartesianProduct(const CartesianExpr<E1, T> &e1,
                                     const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {
#ifndef __CUDACC__
    assert((*this).e1.size() == (*this).e2.size());
#endif
  }

  ANY_DEVICE_INLINE scalar_t operator()(int_t i) const { return e1(i) * e2(i); }

  ANY_DEVICE_INLINE constexpr static int_t size(void) {
    static_assert(E1::size() == E2::size(), "Length mismatch.");
    return E1::size();
  }

private:
  const E1 &e1;
  const E2 &e2;
};

template <class E1, class E2, class T>
class CartesianMinus : public CartesianExpr<CartesianMinus<E1, E2, T>, T> {
public:
  using scalar_t = T;

  ANY_DEVICE_INLINE CartesianMinus(const CartesianExpr<E1, T> &e1,
                                   const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  ANY_DEVICE_INLINE scalar_t operator()(int_t i) const { return e1(i) - e2(i); }

  ANY_DEVICE_INLINE constexpr static int_t size(void) {
    static_assert(E1::size() == E2::size(), "Length mismatch.");
    return E1::size();
  }

private:
  const E1 &e1;
  const E2 &e2;
};

template <class E, class T>
class CartesianScale : public CartesianExpr<CartesianScale<E, T>, T> {
public:
  using scalar_t = T;

  ANY_DEVICE_INLINE CartesianScale(const CartesianExpr<E, T> &e,
                                   const double factor)
      : e(static_cast<const E &>(e)), factor(factor) {}

  ANY_DEVICE_INLINE scalar_t operator()(int_t i) const { return factor * e(i); }

  ANY_DEVICE_INLINE constexpr static int_t size(void) { return E::size(); }

private:
  const E &e;
  const double factor;
};

template <class E1, class E2, class T>
class CartesianDivision
    : public CartesianExpr<CartesianDivision<E1, E2, T>, T> {
public:
  using scalar_t = T;

public:
  ANY_DEVICE_INLINE CartesianDivision(const CartesianExpr<E1, T> &e1,
                                      const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  ANY_DEVICE_INLINE scalar_t operator()(int_t i) const { return e1(i) / e2(i); }

  ANY_DEVICE_INLINE constexpr static int_t size(void) {
    static_assert(E1::size() == E2::size(), "Length mismatch.");
    return E1::size();
  }

private:
  const E1 &e1;
  const E2 &e2;
};

template <class E1, class E2, class T>
class CartesianCross : public CartesianExpr<CartesianCross<E1, E2, T>, T> {
public:
  using scalar_t = T;

  ANY_DEVICE_INLINE CartesianCross(const CartesianExpr<E1, T> &e1,
                                   const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  ANY_DEVICE_INLINE scalar_t operator()(int_t i) const {
    switch (i) {
    case 0:
      return e1(1) * e2(2) - e1(2) * e2(1);
    case 1:
      return e1(2) * e2(0) - e1(0) * e2(2);
    case 2:
      return e1(0) * e2(1) - e1(1) * e2(0);
    default:
      return -42.0;
    }
  }

  ANY_DEVICE_INLINE constexpr static int_t size(void) {
    static_assert(E1::size() == 3, "Cross-product only works for 3-vectors.");
    static_assert(E2::size() == 3, "Cross-product only works for 3-vectors.");

    return E1::size();
  }

private:
  const E1 &e1;
  const E2 &e2;
};

/// Multiply with scalar.
template <class E, class T>
ANY_DEVICE_INLINE const CartesianScale<E, T>
operator*(const double factor, const CartesianExpr<E, T> &e) {
  return CartesianScale<E, T>(e, factor);
}

/// Multiply with scalar.
template <class E, class T>
inline const CartesianScale<E, T> operator*(const CartesianExpr<E, T> &e,
                                            const double factor) {
  return CartesianScale<E, T>(e, factor);
}

/// Divide by scalar
template <class E, class T>
ANY_DEVICE_INLINE const CartesianScale<E, T>
operator/(const CartesianExpr<E, T> &e, const double factor) {
  return CartesianScale<E, T>(e, 1.0 / factor);
}

/// Divide by scalar
template <class E1, class E2, class T>
ANY_DEVICE_INLINE const CartesianDivision<E1, E2, T>
operator/(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianDivision<E1, E2, T>(e1, e2);
}

/// Addition: x + y
template <class E1, class E2, class T>
ANY_DEVICE_INLINE const CartesianPlus<E1, E2, T>
operator+(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianPlus<E1, E2, T>(e1, e2);
}

/// Unary minus: -x
template <class E, class T>
ANY_DEVICE_INLINE const CartesianScale<E, T>
operator-(const CartesianExpr<E, T> &e) {
  return CartesianScale<E, T>(e, -1);
}

/// Subtraction: x - y
template <class E1, class E2, class T>
ANY_DEVICE_INLINE const CartesianMinus<E1, E2, T>
operator-(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianMinus<E1, E2, T>(e1, e2);
}

/// Subtraction: x - y
template <class E1, class E2, class T>
ANY_DEVICE_INLINE const CartesianProduct<E1, E2, T>
operator*(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianProduct<E1, E2, T>(e1, e2);
}

/// Cross (or outer) product.
template <class E1, class E2, class T>
ANY_DEVICE_INLINE const CartesianCross<E1, E2, T>
cross(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianCross<E1, E2, T>(e1, e2);
}

/// Dot (or inner) product.
template <class E1, class E2, class T>
ANY_DEVICE_INLINE double dot(const CartesianExpr<E1, T> &e_,
                             const CartesianExpr<E2, T> &f_) {
  const E1 &e = static_cast<const E1 &>(e_);
  const E2 &f = static_cast<const E2 &>(f_);

  double ef = 0;
  for (int_t i = 0; i < e.size(); ++i) {
    ef += e(i) * f(i);
  }

  return ef;
}

/// Norm induced by the inner product.
template <class T>
ANY_DEVICE_INLINE double norm(const T &e) {
  return std::sqrt(dot(e, e));
}

/// Pointwise almost equal.
/** @param e1  expression to be compared with the second argument
 *  @param e2  expression to be compared with the first argument
 *  @param atol absolute tolerance.
 */
template <class E1, class E2, class T>
ANY_DEVICE_INLINE bool almost_equal(const CartesianExpr<E1, T> &e1,
                                    const CartesianExpr<E2, T> &e2,
                                    double atol) {
  return norm(e2 - e1) < atol;
}

/// Minimum in `e_`.
template <class E, class T>
typename E::scalar_t minimum(const CartesianExpr<E, T> &e_) {
  const E &e = static_cast<const E &>(e_);

  typename E::scalar_t gbl_min = e(0);
  for (int_t i = 1; i < e.size(); ++i) {
    gbl_min = std::min(gbl_min, e(i));
  }

  return gbl_min;
}

/// Maximum in `e_`.
template <class E, class T>
typename E::scalar_t maximum(const CartesianExpr<E, T> &e_) {
  const E &e = static_cast<const E &>(e_);

  typename E::scalar_t gbl_max = e(0);
  for (int_t i = 1; i < e.size(); ++i) {
    gbl_max = std::max(gbl_max, e(i));
  }

  return gbl_max;
}

/// Output
template <class E, class T>
std::ostream &operator<<(std::ostream &os, const CartesianExpr<E, T> &e_) {
  const E &e = static_cast<const E &>(e_);
  os << "[ ";
  for (int_t i = 0; i < e.size(); ++i) {
    os << e(i) << (i == e.size() - 1 ? " ]" : ", ");
  }

  return os;
}

} // namespace zisa

#endif /* end of include guard: CARTESIAN_EXPR_H_KJMKGEYM */
