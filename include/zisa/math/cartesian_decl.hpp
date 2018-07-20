/* Simple fixed size arrays.
 */
#ifndef COORDINATES_DECL_HPP_PCND3ABZ
#define COORDINATES_DECL_HPP_PCND3ABZ

#include "zisa/config.hpp"

#include "zisa/math/cartesian_expr.hpp"
#include "zisa/memory/array_base.hpp"

namespace zisa {
template <int n_vars>
class Cartesian : public CartesianExpr<Cartesian<n_vars>, double> {
public:
  using scalar_t = double;

  /// Empty constructor.
  /** Make CUDA happy.
   */
  ANY_DEVICE_INLINE Cartesian() = default;

  // Copy constructor.
  ANY_DEVICE_INLINE Cartesian(const Cartesian &other) { (*this) = other; }

  // Construct from CRTP expression.
  template <class E>
  explicit Cartesian(const CartesianExpr<E, scalar_t> &e) {
    static_assert(E::size() == size(), "Length mismatch.");
    (*this) = e;
  }

  ANY_DEVICE_INLINE Cartesian(std::initializer_list<double> list) {
    assert(list.size() == n_vars);

    double *dest = (*this).data;
    for (auto &&d : list) {
      *dest = d;
      ++dest;
    }
  }

  ANY_DEVICE_INLINE const scalar_t &operator()(int i) const { return data[i]; }

  ANY_DEVICE_INLINE scalar_t &operator()(int i) { return data[i]; }

  ANY_DEVICE_INLINE const scalar_t &operator[](int i) const { return data[i]; }

  ANY_DEVICE_INLINE scalar_t &operator[](int i) { return data[i]; }

  /// Deep copy.
  ANY_DEVICE_INLINE Cartesian<n_vars> &
  operator=(const Cartesian<n_vars> &other) {
    for (int i = 0; i < n_vars; ++i) {
      data[i] = other.data[i];
    }

    return *this;
  }

  ANY_DEVICE_INLINE Cartesian<n_vars> &
  operator=(std::initializer_list<double> list) {
    assert(list.size() == n_vars);

    double *dest = (*this).data;
    for (auto &&d : list) {
      *dest = d;
      ++dest;
    }

    return *this;
  }

  /// Deep copy.
  template <class E>
  Cartesian<n_vars> &operator=(const CartesianExpr<E, scalar_t> &e_) {
    static_assert(E::size() == size(), "Number of variables don't match!");
    const E &e = static_cast<const E &>(e_);

    for (int i = 0; i < e.size(); ++i) {
      data[i] = e(i);
    }

    return *this;
  }

  /// Add, inplace.
  template <class E>
  void operator+=(const CartesianExpr<E, scalar_t> &e_) {
    const E &e = static_cast<const E &>(e_);

    for (int i = 0; i < e.size(); ++i) {
      data[i] += e(i);
    }
  }

  /// Subtract, inplace.
  template <class E>
  void operator-=(const CartesianExpr<E, scalar_t> &e_) {
    const E &e = static_cast<const E &>(e_);

    for (int i = 0; i < e.size(); ++i) {
      data[i] -= e(i);
    }
  }

  /// Multiply with scalar (inplace)
  inline void operator*=(const double factor) {
    for (int i = 0; i < size(); ++i) {
      data[i] *= factor;
    }
  }

  /// Divide by scalar (inplace)
  inline void operator/=(const double factor) {
    for (int i = 0; i < size(); ++i) {
      data[i] /= factor;
    }
  }

  /// Exactly equal?
  /** See `almost_equal` for equality modulo round-off.
   */
  template <class E>
  bool operator==(const CartesianExpr<E, scalar_t> &e_) const {
    const E &e = static_cast<const E &>(e_);
    for (int i = 0; i < e.size(); ++i) {
      if (data[i] != e(i))
        return false;
    }

    return true;
  }

  /// Not exactly equal.
  template <class E>
  bool operator!=(const CartesianExpr<E, scalar_t> &e) const {
    return !((*this) == e);
  }

  ANY_DEVICE constexpr static int size(void) { return n_vars; }

  ANY_DEVICE static Cartesian<n_vars> zeros(void) {
    Cartesian<n_vars> ret;
    return ret.assign_zero();
  }

  ANY_DEVICE_INLINE Cartesian<n_vars> &assign_zero() {
    return this->assign_constant(0);
  }
  ANY_DEVICE_INLINE Cartesian<n_vars> &assign_constant(const scalar_t &rhs) {
    for (int i = 0; i < n_vars; ++i) {
      (*this)[i] = rhs;
    }

    return *this;
  }

  /// Cartesian orthonormal basis $e_i(j) = delta_{i,j}.$
  ANY_DEVICE static Cartesian<n_vars> unit_vector(int i) {
    Cartesian<n_vars> e = zeros();
    e(i) = 1;

    return e;
  }

protected:
  double data[n_vars];
};

template <int n_vars>
std::ostream &operator<<(std::ostream &os, const Cartesian<n_vars> &x) {

  os << "[ ";
  for (int i = 0; i < x.size(); ++i) {
    os << x[i] << (i == x.size() - 1 ? " ]" : ", ");
  }

  return os;
}

class XY : public Cartesian<2> {
private:
  using super = Cartesian<2>;

public:
  using scalar_t = double;

  using super::Cartesian;
  using super::operator=;
};

template <class E1, class E2, class E3>
double cos_angle(const CartesianExpr<E1, double> &e1_,
                 const CartesianExpr<E2, double> &e2_,
                 const CartesianExpr<E3, double> &e3_) {

  const auto &e1 = static_cast<const E1 &>(e1_);
  const auto &e2 = static_cast<const E2 &>(e2_);
  const auto &e3 = static_cast<const E3 &>(e3_);

  return zisa::dot(normalize(e1 - e2), normalize(e3 - e2));
}

template <class E>
bool isreal(const CartesianExpr<E, double> &e_) {
  const E &e = static_cast<const E &>(e_);

  for (int k = 0; k < e.size(); ++k) {
    if (!isreal(e[k])) {
      return false;
    }
  }

  return true;
}

template <class E>
Cartesian<E::size()> normalize(const CartesianExpr<E, double> &x) {
  return Cartesian<E::size()>(x / norm(x));
}

template <class E>
XY rotate_right(const CartesianExpr<E, double> &e_) {
  const E &e = static_cast<const E &>(e_);
  return {e(1), -e(0)};
}

template <class E>
XY rotate_left(const CartesianExpr<E, double> &e_) {
  const E &e = static_cast<const E &>(e_);
  return {-e(1), e(0)};
}

/// Save an array of `XY`.
template <class Indexing, class Array, class Shape>
static void save(HDF5Writer &writer,
                 const array_base<XY, Indexing, Array, Shape> &arr,
                 const std::string &tag) {
  XY const *const data = arr.raw();
  const auto &dims = arr.shape;

  HDF5DataType data_type = make_hdf5_data_type<double>();

  constexpr int rank = Shape::size() + 1;
  hsize_t h5_dims[rank];
  for (int i = 0; i < rank - 1; ++i) {
    h5_dims[i] = hsize_t(dims(i)); // size of (i, j, k) axes
  }
  h5_dims[rank - 1] = 2; // number of space components

  writer.write_array((double *)data, data_type, tag, rank, h5_dims);
}

// /// Save an array of `XY`.
// template <class LinearIndex>
// void load_impl(HDF5Reader &reader,
//                ArrayBase<XY, LinearIndex> &array,
//                const std::string &tag) {
//   constexpr int n_dims = LinearIndex::n_dims;
//   hsize_t dims[n_dims + 1];

//   HDF5DataType data_type = make_hdf5_data_type<double>();
//   XY *raw_data
//       = (XY *)reader.read_array(data_type, tag, n_dims + 1, dims);

//   Shape<n_dims> shape;
//   for (int i = 0; i < n_dims; ++i) {
//     shape(i) = int(dims[i]);
//   }

//   array = Array<XY, n_dims>(raw_data, shape);
// }

} // namespace zisa

#endif /* end of include guard: COORDINATES_H_PCND3ABZ */
