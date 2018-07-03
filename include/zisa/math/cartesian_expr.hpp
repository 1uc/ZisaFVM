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

  __host__ __device__ __inline__ CartesianPlus(const CartesianExpr<E1, T> &e1,
                                               const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  __host__ __device__ __inline__ scalar_t operator()(int i) const {
    return e1(i) + e2(i);
  }

  __host__ __device__ __inline__ constexpr static int size(void) {
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

  __host__
      __device__ __inline__ CartesianProduct(const CartesianExpr<E1, T> &e1,
                                             const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {
#ifndef __CUDACC__
    assert((*this).e1.size() == (*this).e2.size());
#endif
  }

  __host__ __device__ __inline__ scalar_t operator()(int i) const {
    return e1(i) * e2(i);
  }

  __host__ __device__ __inline__ constexpr static int size(void) {
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

  __host__ __device__ __inline__ CartesianMinus(const CartesianExpr<E1, T> &e1,
                                                const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  __host__ __device__ __inline__ scalar_t operator()(int i) const {
    return e1(i) - e2(i);
  }

  __host__ __device__ __inline__ constexpr static int size(void) {
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

  __host__ __device__ __inline__ CartesianScale(const CartesianExpr<E, T> &e,
                                                const double factor)
      : e(static_cast<const E &>(e)), factor(factor) {}

  __host__ __device__ __inline__ scalar_t operator()(int i) const {
    return factor * e(i);
  }

  __host__ __device__ __inline__ constexpr static int size(void) {
    return E::size();
  }

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
  __host__
      __device__ __inline__ CartesianDivision(const CartesianExpr<E1, T> &e1,
                                              const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  __host__ __device__ __inline__ scalar_t operator()(int i) const {
    return e1(i) / e2(i);
  }

  __host__ __device__ __inline__ constexpr static int size(void) {
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

  __host__ __device__ __inline__ CartesianCross(const CartesianExpr<E1, T> &e1,
                                                const CartesianExpr<E2, T> &e2)
      : e1(static_cast<const E1 &>(e1)), e2(static_cast<const E2 &>(e2)) {}

  __host__ __device__ __inline__ scalar_t operator()(int i) const {
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

  __host__ __device__ __inline__ constexpr static int size(void) {
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
__host__ __device__ __inline__ const CartesianScale<E, T>
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
__host__ __device__ __inline__ const CartesianScale<E, T>
operator/(const CartesianExpr<E, T> &e, const double factor) {
  return CartesianScale<E, T>(e, 1.0 / factor);
}

/// Divide by scalar
template <class E1, class E2, class T>
__host__ __device__ __inline__ const CartesianDivision<E1, E2, T>
operator/(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianDivision<E1, E2, T>(e1, e2);
}

/// Addition: x + y
template <class E1, class E2, class T>
__host__ __device__ __inline__ const CartesianPlus<E1, E2, T>
operator+(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianPlus<E1, E2, T>(e1, e2);
}

/// Unary minus: -x
template <class E, class T>
__host__ __device__ __inline__ const CartesianScale<E, T>
operator-(const CartesianExpr<E, T> &e) {
  return CartesianScale<E, T>(e, -1);
}

/// Subtraction: x - y
template <class E1, class E2, class T>
__host__ __device__ __inline__ const CartesianMinus<E1, E2, T>
operator-(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianMinus<E1, E2, T>(e1, e2);
}

/// Subtraction: x - y
template <class E1, class E2, class T>
__host__ __device__ __inline__ const CartesianProduct<E1, E2, T>
operator*(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianProduct<E1, E2, T>(e1, e2);
}

/// Cross (or outer) product.
template <class E1, class E2, class T>
__host__ __device__ __inline__ const CartesianCross<E1, E2, T>
cross(const CartesianExpr<E1, T> &e1, const CartesianExpr<E2, T> &e2) {
  return CartesianCross<E1, E2, T>(e1, e2);
}

/// Dot (or inner) product.
template <class E1, class E2, class T>
__host__ __device__ __inline__ double dot(const CartesianExpr<E1, T> &e_,
                                          const CartesianExpr<E2, T> &f_) {
  const E1 &e = static_cast<const E1 &>(e_);
  const E2 &f = static_cast<const E2 &>(f_);

  double ef = 0;
  for (int i = 0; i < e.size(); ++i) {
    ef += e(i) * f(i);
  }

  return ef;
}

/// Norm induced by the inner product.
template <class T>
__host__ __device__ __inline__ double norm(const T &e) {
  return std::sqrt(dot(e, e));
}

/// Pointwise almost equal.
/** @param e1  expression to be compared with the second argument
 *  @param e2  expression to be compared with the first argument
 *  @param atol absolute tolerance.
 */
template <class E1, class E2, class T>
__host__ __device__ __inline__ bool almost_equal(const CartesianExpr<E1, T> &e1,
                                                 const CartesianExpr<E2, T> &e2,
                                                 double atol) {
  return norm(e2 - e1) < atol;
}

/// Minimum in `e_`.
template <class E, class T>
typename E::scalar_t minimum(const CartesianExpr<E, T> &e_) {
  const E &e = static_cast<const E &>(e_);

  typename E::scalar_t gbl_min = e(0);
  for (int i = 1; i < e.size(); ++i) {
    gbl_min = std::min(gbl_min, e(i));
  }

  return gbl_min;
}

/// Maximum in `e_`.
template <class E, class T>
typename E::scalar_t maximum(const CartesianExpr<E, T> &e_) {
  const E &e = static_cast<const E &>(e_);

  typename E::scalar_t gbl_max = e(0);
  for (int i = 1; i < e.size(); ++i) {
    gbl_max = std::max(gbl_max, e(i));
  }

  return gbl_max;
}

/// Output
template <class E, class T>
std::ostream &operator<<(std::ostream &os, const CartesianExpr<E, T> &e_) {
  const E &e = static_cast<const E &>(e_);
  os << "[ ";
  for (int i = 0; i < e.size(); ++i) {
    os << e(i) << (i == e.size() - 1 ? " ]" : ", ");
  }

  return os;
}

} // namespace zisa

#endif /* end of include guard: CARTESIAN_EXPR_H_KJMKGEYM */
