#ifndef GRAVITY_H_WB0OY
#define GRAVITY_H_WB0OY

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

template <class GravityBase, class Alignment>
class Gravity {
public:
  Gravity() = default;
  Gravity(const GravityBase &gravity, const Alignment &alignment)
      : gravity(gravity), alignment(alignment) {}

  ANY_DEVICE_INLINE double phi(const XY &x) const {
    double chi = alignment.coordinate(x);
    return gravity.phi(chi);
  }

  ANY_DEVICE_INLINE double dphi_dx(const XY &x, int dir) const {
    double chi = alignment.coordinate(x);
    return gravity.dphi_dx(chi) * alignment.dx(x, dir);
  }

  ANY_DEVICE_INLINE XY grad_phi(const XY &x) const {
    return {dphi_dx(x, 0), dphi_dx(x, 1)};
  }

  ANY_DEVICE_INLINE double norm_grad_phi(const XY &x) const {
    return zisa::norm(grad_phi(x));
  }

protected:
  GravityBase gravity;
  Alignment alignment;
};

class RadialAlignment {
public:
  ANY_DEVICE_INLINE double coordinate(const XY &xy) const {
    return zisa::norm(xy);
  }

  ANY_DEVICE_INLINE double dx(const XY &xy, int dir) const {
    double r = zisa::norm(xy);
    return xy[dir] / (r + epsilon);
  }

private:
  double epsilon = 1e-50;
};

class AxialAlignment {
public:
  ANY_DEVICE_INLINE AxialAlignment() = default;
  ANY_DEVICE_INLINE AxialAlignment(const XY &axis) : axis(axis){};

  ANY_DEVICE_INLINE double coordinate(const XY &xy) const {
    return zisa::dot(xy, axis);
  }
  ANY_DEVICE_INLINE double dx(const XY &, int dir) const { return axis[dir]; }

private:
  XY axis;
};

class ConstantGravity {
public:
  ANY_DEVICE_INLINE ConstantGravity() = default;
  ANY_DEVICE_INLINE ConstantGravity(double gravity) : gravity(gravity) {}

  ANY_DEVICE_INLINE double phi(double chi) const;
  ANY_DEVICE_INLINE double dphi_dx(double chi) const;

private:
  double gravity;
};

class ConstantGravityRadial : public Gravity<ConstantGravity, RadialAlignment> {
private:
  using super = Gravity<ConstantGravity, RadialAlignment>;

public:
  ConstantGravityRadial() = default;
  ConstantGravityRadial(double gravity) : super(gravity, RadialAlignment()){};
};

class ConstantGravityAxial : public Gravity<ConstantGravity, AxialAlignment> {
private:
  using super = Gravity<ConstantGravity, AxialAlignment>;

public:
  using super::super;
};

class PointMassGravity {
public:
  PointMassGravity() = default;
  PointMassGravity(double G, double M, double X);

  ANY_DEVICE_INLINE double phi(double chi) const;
  ANY_DEVICE_INLINE double dphi_dx(double chi) const;

private:
  double GM;
  double X;
};

class PointMassGravityRadial
    : public Gravity<PointMassGravity, RadialAlignment> {
private:
  using super = Gravity<PointMassGravity, RadialAlignment>;

public:
  PointMassGravityRadial() = default;
  PointMassGravityRadial(double G, double M, double X)
      : super({G, M, X}, RadialAlignment()){};
};

class PointMassGravityAxial : public Gravity<PointMassGravity, AxialAlignment> {
private:
  using super = Gravity<PointMassGravity, AxialAlignment>;

public:
  PointMassGravityAxial() = default;
  using super::super;
};

class PointMassGravityXAligned : public PointMassGravityAxial {
private:
  using super = PointMassGravityAxial;

public:
  PointMassGravityXAligned() = default;
  PointMassGravityXAligned(double G, double M, double X)
      : super({G, M, X}, XY{1.0, 0.0}){};
};

class PointMassGravityYAligned : public PointMassGravityAxial {
private:
  using super = PointMassGravityAxial;

public:
  PointMassGravityYAligned() = default;
  PointMassGravityYAligned(double G, double M, double X)
      : super({G, M, X}, XY{0.0, 1.0}){};
};

class PolytropeGravity {
public:
  PolytropeGravity() = default;
  __host__ PolytropeGravity(double rhoC, double K, double G);
  ANY_DEVICE_INLINE double phi(double chi) const;
  ANY_DEVICE_INLINE double dphi_dx(double chi) const;

  ANY_DEVICE_INLINE double alpha() const;

private:
  double rhoC = 1.0;
  double K = 1.0;
  double G = 1.0;
  double eps = 1e-50;
};

class PolytropeGravityRadial
    : public Gravity<PolytropeGravity, RadialAlignment> {
private:
  using super = Gravity<PolytropeGravity, RadialAlignment>;

public:
  PolytropeGravityRadial() = default;
  PolytropeGravityRadial(double rhoC, double K, double G)
      : super({rhoC, K, G}, RadialAlignment{}) {}

  ANY_DEVICE_INLINE double alpha() const;
};

class SphericalGravity {
private:
  class Poly1D {
  public:
    Poly1D() = default;
    Poly1D(double a1, double a2, double a3, double b)
        : a1(a1), a2(a2), a3(a3), b(b) {}

    ANY_DEVICE_INLINE double operator()(double r) const {
      return ((a3 * r + a2) * r + a1) * r + b * (r - 1);
    }

    ANY_DEVICE_INLINE double dr(double r) const {
      return (3 * a3 * r + 2 * a2) * r + a1 + b;
    }

    ANY_DEVICE_INLINE double drr(double r) const { return 6 * a3 * r + 2 * a2; }

  private:
    double a1;
    double a2;
    double a3;
    double b;
  };

public:
  SphericalGravity() = default;
  SphericalGravity(const std::vector<double> &domain,
                   const std::vector<double> &phi);

  ANY_DEVICE_INLINE double phi(double r) const;
  ANY_DEVICE_INLINE double dphi_dx(double r) const;

private:
  ANY_DEVICE_INLINE int index(double r) const;
  ANY_DEVICE_INLINE double radii(int i) const;

  ANY_DEVICE_INLINE Poly1D make_poly(double f_i,
                                     double df_i,
                                     double ddf_i,
                                     double f_ip1) const;

private:
  std::vector<double> domain;
  int n_cells;
  double dr;

  std::vector<double> phi_points;
};

class RadialGravity : public Gravity<SphericalGravity, RadialAlignment> {
private:
  using super = Gravity<SphericalGravity, RadialAlignment>;

public:
  RadialGravity() = default;
  RadialGravity(const std::vector<double> &domain,
                const std::vector<double> &phi)
      : super({domain, phi}, RadialAlignment{}) {}
};

class NoGravity {
public:
  template <class... Args>
  double phi(const Args &... /* args */) const {
    return 0.0;
  }

  template <class... Args>
  double dphi_dx(const Args &... /* args */) const {
    return 0.0;
  }
};

} // namespace zisa

#endif /* end of include guard */
