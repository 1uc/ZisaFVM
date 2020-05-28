#ifndef CONE_H_CKYAW
#define CONE_H_CKYAW

#include <zisa/config.hpp>

#include <vector>
#include <zisa/math/cartesian.hpp>

namespace zisa {

class Region {
public:
  virtual ~Region() = default;
  virtual bool is_inside(const XYZ &x) const = 0;
};

/// Two dimensional wedge.
class TriangularCone : public Region {
public:
  /// The point of the wedge is A; ABC is right-handed.
  TriangularCone(const XYZ &A, const XYZ &B, const XYZ &C);

  bool is_inside(const XYZ &x) const override;

private:
  XYZ A;
  XYZ dB;
  XYZ dC;
};

/// An unbounded tetrahedron.
class TetrahedralCone : public Region {
public:
  /// The point of the 3D wedge is A; BCD points away from A.
  TetrahedralCone(const XYZ &A, const XYZ &B, const XYZ &C, const XYZ &D);

  bool is_inside(const XYZ &x) const override;

private:
  XYZ A;
  XYZ dB;
  XYZ dC;
  XYZ dD;
};

/// A ball of finite radius.
class Ball : public Region {
public:
  Ball(const XYZ &x_center, double radius);

  bool is_inside(const XYZ &x) const override;

private:
  XYZ x_center;
  double radius;
};

/// Always true.
/** Useful if the region is a cone with angle 2 pi, e.g. everything.
 */
class FullSphere : public Region {
public:
  bool is_inside(const XYZ &x) const override;
};

class UnionOfRegions : public Region {
public:
  bool is_inside(const XYZ &x) const override;

  void add(std::shared_ptr<Region> region);

private:
  std::vector<std::shared_ptr<Region>> regions;
};

} // namespace zisa
#endif /* end of include guard */
