#ifndef GRAVITY_H_WB0OY
#define GRAVITY_H_WB0OY

#include <limits>

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class GravityBase, class Alignment>
class Gravity {
public:
  Gravity() = default;
  Gravity(GravityBase gravity, Alignment alignment)
      : gravity(std::move(gravity)), alignment(std::move(alignment)) {}

  ANY_DEVICE_INLINE double phi(const XYZ &x) const {
    double chi = alignment.coordinate(x);
    return gravity.phi(chi);
  }

  ANY_DEVICE_INLINE double dphi_dx(const XYZ &x, int_t dir) const {
    double chi = alignment.coordinate(x);
    return gravity.dphi_dx(chi) * alignment.dx(x, dir);
  }

  ANY_DEVICE_INLINE XYZ grad_phi(const XYZ &x) const {
    return {dphi_dx(x, 0), dphi_dx(x, 1), dphi_dx(x, 2)};
  }

  ANY_DEVICE_INLINE double norm_grad_phi(const XYZ &x) const {
    return zisa::norm(grad_phi(x));
  }

  std::string str() const {
    return string_format("%s with %s",
                         type_name<GravityBase>().c_str(),
                         type_name<Alignment>().c_str());
  }

  template <class G, class A>
  friend void save(HDF5Writer &writer, const Gravity<G, A> &g);

protected:
  GravityBase gravity;
  Alignment alignment;
};

template <class GravityBase, class Alignment>
void save(HDF5Writer &writer, const Gravity<GravityBase, Alignment> &gravity) {
  writer.open_group("gravity");
  save(writer, gravity.gravity);
  writer.switch_group("alignment");
  save(writer, gravity.alignment);
  writer.close_group();
}

class RadialAlignment {
public:
  ANY_DEVICE_INLINE double coordinate(const XYZ &xy) const {
    return zisa::norm(xy);
  }

  ANY_DEVICE_INLINE double dx(const XYZ &xy, int_t dir) const {
    double r = zisa::norm(xy);
    return xy[dir] / (r + epsilon);
  }

  friend void save(HDF5Writer &writer, const RadialAlignment &alignment);

private:
  double epsilon = 1e-50;
};

void save(HDF5Writer &writer, const RadialAlignment &alignment);

class AxialAlignment {
public:
  ANY_DEVICE_INLINE AxialAlignment() = default;
  ANY_DEVICE_INLINE AxialAlignment(const XYZ &axis) : axis(axis){};

  ANY_DEVICE_INLINE double coordinate(const XYZ &xy) const {
    return zisa::dot(xy, axis);
  }
  ANY_DEVICE_INLINE double dx(const XYZ &, int_t dir) const {
    return axis[dir];
  }

  friend void save(HDF5Writer &writer, const AxialAlignment &alignment);

private:
  XYZ axis;
};

void save(HDF5Writer &writer, const AxialAlignment &alignment);

class ConstantGravity {
public:
  ANY_DEVICE_INLINE ConstantGravity() = default;
  ANY_DEVICE_INLINE ConstantGravity(double gravity) : gravity(gravity) {}

  ANY_DEVICE_INLINE double phi(double chi) const;
  ANY_DEVICE_INLINE double dphi_dx(double chi) const;

  friend void save(HDF5Writer &writer, const ConstantGravity &gravity);

private:
  double gravity;
};

void save(HDF5Writer &writer, const ConstantGravity &gravity);

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

  friend void save(HDF5Writer &writer, const PointMassGravity &alignment);

private:
  double GM;
  double X;
};

void save(HDF5Writer &writer, const PointMassGravity &gravity);

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
      : super({G, M, X}, XYZ::unit_vector(0)){};
};

class PointMassGravityYAligned : public PointMassGravityAxial {
private:
  using super = PointMassGravityAxial;

public:
  PointMassGravityYAligned() = default;
  PointMassGravityYAligned(double G, double M, double X)
      : super({G, M, X}, XYZ::unit_vector(1)){};
};

class PolytropeGravity {
public:
  PolytropeGravity() = default;
  inline PolytropeGravity(double rhoC, double K, double G);

  ANY_DEVICE_INLINE double phi(double chi) const;
  ANY_DEVICE_INLINE double dphi_dx(double chi) const;

  ANY_DEVICE_INLINE RhoEntropy rhoK_center() const;
  ANY_DEVICE_INLINE double alpha(double chi) const;

  friend void save(HDF5Writer &writer, const PolytropeGravity &gravity);

private:
  double rhoC = 1.0;
  double K = 1.0;
  double G = 1.0;
  double eps = std::numeric_limits<double>::min();
};

class PolytropeGravityRadial
    : public Gravity<PolytropeGravity, RadialAlignment> {
private:
  using super = Gravity<PolytropeGravity, RadialAlignment>;

public:
  PolytropeGravityRadial() = default;
  PolytropeGravityRadial(double rhoC, double K, double G)
      : super(PolytropeGravity(rhoC, K, G), RadialAlignment{}) {}

  ANY_DEVICE_INLINE double alpha(double chi) const;
  ANY_DEVICE_INLINE RhoEntropy rhoK_center() const;
};

class PolytropeGravityWithJump {
public:
  PolytropeGravityWithJump() = default;
  inline PolytropeGravityWithJump(
      double r_crit, double rhoC, double K_inner, double K_outer, double G);
  ANY_DEVICE_INLINE double phi(double chi) const;
  ANY_DEVICE_INLINE double dphi_dx(double chi) const;

  ANY_DEVICE_INLINE double alpha(double chi) const;

  friend void save(HDF5Writer &writer, const PolytropeGravityWithJump &gravity);

private:
  double r_crit;
  PolytropeGravity inner;
  PolytropeGravity outer;
};

void save(HDF5Writer &writer, const PolytropeGravityWithJump &gravity);

class PolytropeGravityWithJumpRadial
    : public Gravity<PolytropeGravityWithJump, RadialAlignment> {
private:
  using super = Gravity<PolytropeGravityWithJump, RadialAlignment>;

public:
  PolytropeGravityWithJumpRadial() = default;

  PolytropeGravityWithJumpRadial(
      double r_crit, double rhoC, double K_inner, double K_outer, double G)
      : super({r_crit, rhoC, K_inner, K_outer, G}, RadialAlignment{}) {}

  ANY_DEVICE_INLINE double alpha(double chi) const;
};

class SphericalGravity {
private:
public:
  SphericalGravity() = default;
  SphericalGravity(array<double, 1> radii, array<double, 1> phi);

  inline double phi(double r) const;
  inline double dphi_dx(double r) const;

private:
  inline int_t index(double r) const;
  inline double radii(int_t i) const;

  friend void save(HDF5Writer &writer, const SphericalGravity &gravity);

private:
  std::shared_ptr<array<double, 1>> radii_;
  std::shared_ptr<array<double, 1>> phi_;
  int_t n_cells;
};

void save(HDF5Writer &writer, const SphericalGravity &gravity);

class RadialGravity : public Gravity<SphericalGravity, RadialAlignment> {
private:
  using super = Gravity<SphericalGravity, RadialAlignment>;

public:
  RadialGravity() = default;
  RadialGravity(array<double, 1> radii, array<double, 1> phi)
      : super({std::move(radii), std::move(phi)}, RadialAlignment{}) {}
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

  inline std::string str() const { return "No gravity."; }
};

} // namespace zisa

#endif /* end of include guard */
