// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef POLY2D_IMPL_H_PP21T
#define POLY2D_IMPL_H_PP21T

#include "poly2d_decl.hpp"

namespace zisa {

template <int NDIMS>
constexpr ANY_DEVICE_INLINE int_t poly_dof(int deg) {
  if constexpr (NDIMS == 2) {
    return int_t(((deg + 1) * (deg + 2)) / 2);
  } else {
    return int_t(((deg + 1) * (deg + 2) * (deg + 3)) / 6);
  }
}

ANY_DEVICE int_t poly_dof(int deg, int n_dims);

template <int NDIMS>
int poly_degree(int_t n_coeffs) {

  LOG_ERR_IF(n_coeffs == 0, "Too few coefficients.");

  for (int d = 1; true; ++d) {
    if (poly_dof<NDIMS>(d) > n_coeffs) {
      return d - 1;
    }
  }
}

constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b) {
  auto n = a + b;
  return int_t(((n + 1) * n) / 2 + int_t(b));
}

constexpr ANY_DEVICE_INLINE int_t poly_index(int a, int b, int c) {
  return poly_dof<3>(a + b + c - 1) + poly_index(b, c);
}

template <int NCOEFFS, int NVARS>
int PolyND<NCOEFFS, NVARS>::max_degree() const {
  return poly_degree(n_coeffs(), n_dims());
}

template <int NCOEFFS, int NVARS>
int PolyND<NCOEFFS, NVARS>::degree() const {
  return degree_;
}

template <int NCOEFFS, int NVARS>
int_t PolyND<NCOEFFS, NVARS>::dof() const {
  return poly_dof(degree(), n_dims());
}

template <int NCOEFFS, int NVARS>
constexpr int_t PolyND<NCOEFFS, NVARS>::n_coeffs() {
  return NCOEFFS;
}

template <int NCOEFFS, int NVARS>
constexpr int_t PolyND<NCOEFFS, NVARS>::n_vars() {
  return NVARS;
}

template <int NCOEFFS, int NVARS>
int PolyND<NCOEFFS, NVARS>::n_dims() const {
  return n_dims_;
}

template <int NCOEFFS, int NVARS>
PolyND<NCOEFFS, NVARS>::PolyND() : degree_(0), n_dims_(0) {
  std::fill(coeffs, coeffs + n_coeffs() * n_vars(), 0.0);
  std::fill(moments, moments + n_coeffs(), 0.0);
  x_center_ = XYZ::zeros();
  reference_length_ = 1.0;
}

template <int NCOEFFS, int NVARS>
PolyND<NCOEFFS, NVARS>::PolyND(int degree,
                               const array<double, 1> &moments,
                               const XYZ &x_center,
                               double reference_length,
                               int n_dims)
    : degree_(degree),
      n_dims_(n_dims),
      x_center_(x_center),
      reference_length_(reference_length) {

  std::fill(this->coeffs, this->coeffs + n_coeffs() * n_vars(), 0.0);

  static_assert(n_coeffs() >= 3, "Unusually low degree; fix code below.");
  this->moments[0] = this->moments[1] = this->moments[2] = 0.0;

  if (std::distance(moments.begin(), moments.end()) > 3) {
    std::copy(moments.begin() + 3, moments.end(), this->moments + 3);
  }
}

template <int NCOEFFS, int NVARS>
PolyND<NCOEFFS, NVARS>::PolyND(
    int degree,
    const std::initializer_list<double> &moments_list,
    const XYZ &x_center,
    double reference_length,
    int n_dims)
    : degree_(degree),
      n_dims_(n_dims),
      x_center_(x_center),
      reference_length_(reference_length) {

  std::fill(this->coeffs, this->coeffs + n_coeffs() * n_vars(), 0.0);

  static_assert(n_coeffs() >= 3, "Unusually low degree; fix code below.");
  this->moments[0] = this->moments[1] = this->moments[2] = 0.0;

  if (std::distance(moments_list.begin(), moments_list.end()) > 3) {
    std::copy(moments_list.begin() + 3, moments_list.end(), moments + 3);
  }
}

template <int NCOEFFS, int NVARS>
PolyND<NCOEFFS, NVARS>::PolyND(
    const std::initializer_list<double> &coeffs_list,
    const std::initializer_list<double> &moments_list,
    const XYZ &x_center,
    double reference_length,
    int n_dims)
    : degree_(poly_degree(coeffs_list.size(), n_dims)),
      n_dims_(n_dims),
      x_center_(x_center),
      reference_length_(reference_length) {

  assert(n_vars() == 1);
  assert(coeffs_list.size() <= n_coeffs());
  assert(coeffs_list.size() == moments_list.size());
  assert(poly_dof(degree(), n_dims) == coeffs_list.size());

  std::fill(coeffs, coeffs + n_vars() * n_coeffs(), 0.0);
  std::copy(coeffs_list.begin(), coeffs_list.end(), coeffs);

  std::fill(moments, moments + n_coeffs(), 0.0);
  std::copy(moments_list.begin(), moments_list.end(), moments);
}

template <int NCOEFFS, int NVARS>
template <class E>
PolyND<NCOEFFS, NVARS>::PolyND(const PolynomialCRTP<E> &e_) {
  const auto &e = static_cast<const E &>(e_);

  static_assert(E::n_coeffs() == NCOEFFS, "Number of coefficients mismatch.");
  static_assert(E::n_vars() == NVARS, "Number of variables mismatch.");
  deep_copy(e);
}

template <int NCOEFFS, int NVARS>
double *PolyND<NCOEFFS, NVARS>::coeffs_ptr() {
  return coeffs;
}

template <int NCOEFFS, int NVARS>
template <class E>
void PolyND<NCOEFFS, NVARS>::deep_copy(const PolynomialCRTP<E> &e_) {
  const E &e = static_cast<const E &>(e_);

  static_assert(E::n_coeffs() == NCOEFFS, "Number of coefficients mismatch.");
  static_assert(E::n_vars() == NVARS, "Number of variables mismatch.");

  degree_ = e.degree();

  constexpr int_t n = n_coeffs() * n_vars();
  for (int_t i = 0; i < n; ++i) {
    coeffs[i] = e.a(i);
  }

  for (int_t i = 0; i < n_coeffs(); ++i) {
    moments[i] = e.c(i);
  }

  x_center_ = e.x_center();
  reference_length_ = e.reference_length();
  n_dims_ = e.n_dims();
}

template <int NCOEFFS, int NVARS>
double PolyND<NCOEFFS, NVARS>::a(int_t i) const {
  return coeffs[i];
}

template <int NCOEFFS, int NVARS>
double &PolyND<NCOEFFS, NVARS>::a(int_t i) {
  return coeffs[i];
}

template <int NCOEFFS, int NVARS>
double PolyND<NCOEFFS, NVARS>::c(int_t i) const {
  return moments[i];
}

template <int NCOEFFS, int NVARS>
const XYZ &PolyND<NCOEFFS, NVARS>::x_center() const {
  return x_center_;
}

template <int NCOEFFS, int NVARS>
double PolyND<NCOEFFS, NVARS>::reference_length() const {
  return reference_length_;
}

template <int NCOEFFS, int NVARS>
template <class E>
void PolyND<NCOEFFS, NVARS>::operator=(const PolynomialCRTP<E> &e_) {
  static_assert(E::n_coeffs() == NCOEFFS, "Number of coefficients mismatch.");
  static_assert(E::n_vars() == NVARS, "Number of variables mismatch.");

  const E &e = static_cast<const E &>(e_);
  this->deep_copy(e);
}

template <int NCOEFFS, int NVARS>
template <class E>
void PolyND<NCOEFFS, NVARS>::operator+=(const PolynomialCRTP<E> &e_) {
  static_assert(E::n_coeffs() == NCOEFFS, "Number of coefficients mismatch.");
  static_assert(E::n_vars() == NVARS, "Number of variables mismatch.");

  const E &e = static_cast<const E &>(e_);
  *this = *this + e;
}

template <int NCOEFFS, int NVARS>
template <class E>
void PolyND<NCOEFFS, NVARS>::operator-=(const PolynomialCRTP<E> &e_) {
  static_assert(E::n_coeffs() == NCOEFFS, "Number of coefficients mismatch.");
  static_assert(E::n_vars() == NVARS, "Number of variables mismatch.");

  const E &e = static_cast<const E &>(e_);
  *this = *this - e;
}

template <int NCOEFFS, int NVARS>
void PolyND<NCOEFFS, NVARS>::operator*=(double alpha) {
  *this = alpha * (*this);
}

template <int NCOEFFS, int NVARS>
void PolyND<NCOEFFS, NVARS>::operator/=(double alpha) {
  *this *= 1.0 / alpha;
}

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> PolyND<MAX_DEGREE, NVARS>::operator()(const XYZ &xyz) const {
  if (n_dims() == 2) {
    return eval_2d(xyz);
  } else {
    return eval_3d(xyz);
  }
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t PolyND<MAX_DEGREE, NVARS>::idx(int k, int l) {
  return poly_index(k, l);
}

template <int MAX_DEGREE, int NVARS>
constexpr int_t PolyND<MAX_DEGREE, NVARS>::idx(int k, int l, int m) {
  return poly_index(k, l, m);
}

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> PolyND<MAX_DEGREE, NVARS>::eval_2d(const XYZ &xy) const {
  auto [x, y, _] = XYZ((xy - this->x_center()) / this->reference_length());

  auto d = this->degree();

  Cartesian<NVARS> px(0.0);
  double pow_x_kx = 1.0;
  for (int kx = 0; kx <= d; ++kx) {
    double pow_y_ky = 1.0;
    for (int ky = 0; ky <= d - kx; ++ky) {
      int_t i = idx(kx, ky);
      for (int_t k = 0; k < NVARS; ++k) {
        px[k] += this->a(i * NVARS + k) * (pow_x_kx * pow_y_ky - this->c(i));
      }
      pow_y_ky *= y;
    }

    pow_x_kx *= x;
  }

  return px;
}

template <int MAX_DEGREE, int NVARS>
Cartesian<NVARS> PolyND<MAX_DEGREE, NVARS>::eval_3d(const XYZ &xyz) const {
  auto [x, y, z] = XYZ((xyz - this->x_center()) / this->reference_length());

  auto d = this->degree();

  Cartesian<NVARS> px(0.0);
  double pow_x_kx = 1.0;
  for (int kx = 0; kx <= d; ++kx) {
    double pow_y_ky = 1.0;
    for (int ky = 0; ky <= d - kx; ++ky) {
      double pow_z_kz = 1.0;

      for (int kz = 0; kz <= d - kx - ky; ++kz) {
        int_t i = idx(kx, ky, kz);
        for (int_t k = 0; k < NVARS; ++k) {
          double ak = this->a(i * NVARS + k);
          px[k] += ak * (pow_x_kx * pow_y_ky * pow_z_kz - this->c(i));
        }

        pow_z_kz *= z;
      }
      pow_y_ky *= y;
    }
    pow_x_kx *= x;
  }

  return px;
}

template <int NCOEFFS, int NVARS>
Cartesian<NVARS> smoothness_indicator(const PolyND<NCOEFFS, NVARS> &p) {
  auto beta = Cartesian<NVARS>(0.0);

  int_t n = p.dof();
  for (int_t i = poly_dof(0, p.n_dims()); i < n; ++i) {
    for (int_t k = 0; k < p.n_vars(); ++k) {
      beta[k] += zisa::pow<2>(p.a(i * NVARS + k));
    }
  }

  return beta;
}

template <int NCOEFFS, int NVARS>
std::ostream &operator<<(std::ostream &os, const PolyND<NCOEFFS, NVARS> &poly) {
  os << format_as_list(poly.coeffs, poly.dof()) << "\n";
  os << format_as_list(poly.moments, poly.dof()) << "\n";

  return os;
}

template <int_t NDIMS>
class PolyIndexRange {};

template <>
class PolyIndexRange<2> {
public:
  class EndIterator {
  public:
    EndIterator(int deg) : deg(deg) {}

  public:
    int deg;
  };

  class Iterator {
  public:
    std::array<int, 2> operator*() const { return {d - j, j}; }

    void operator++() {
      ++j;

      if (j > d) {
        ++d;
        j = 0;
      }
    }

    int level() const { return d; }

    bool operator==(const EndIterator &end_it) const {
      return level() > end_it.deg;
    }

    bool operator!=(const EndIterator &end_it) const {
      return !((*this) == end_it);
    }

  private:
    int d = 0;
    int j = 0;
  };

  PolyIndexRange(int max_deg) : max_deg(max_deg) {}

  Iterator begin() const { return Iterator(); }
  EndIterator end() const { return EndIterator(max_deg); }

private:
  int max_deg;
};

template <>
class PolyIndexRange<3> {
private:
  static constexpr int n_dims() { return 3; }

  class EndIterator {
  public:
    EndIterator(int deg) : deg(deg) {}

  public:
    int deg;
  };

  class Iterator {
  public:
    void operator++() {
      ++it_2d;

      if (it_2d.level() > d) {
        ++d;
        it_2d = PolyIndexRange<2>::Iterator();
      }
    }

    std::array<int, 3> operator*() const {
      auto ij = *it_2d;
      return {d - it_2d.level(), ij[0], ij[1]};
    }

    int level() const { return d; }

    bool operator==(const EndIterator &end_it) const {
      return level() > end_it.deg;
    }

    bool operator!=(const EndIterator &end_it) const {
      return !((*this) == end_it);
    }

  private:
    int d = 0;
    PolyIndexRange<2>::Iterator it_2d;
  };

public:
  PolyIndexRange(int max_deg) : max_deg(max_deg) {}

  Iterator begin() const { return Iterator(); }
  EndIterator end() const { return EndIterator(max_deg); }

public:
  int max_deg;
};

} // namespace zisa
#endif /* end of include guard */
