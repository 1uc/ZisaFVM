#include <zisa/math/cone.hpp>
namespace zisa {

TriangularCone::TriangularCone(const XYZ &A, const XYZ &B, const XYZ &C)
    : A(A), dB(XYZ(B - A)), dC(XYZ(C - A)) {}

bool TriangularCone::is_inside(const XYZ &x) const {
  auto dx = XYZ(x - A);
  return zisa::cross(dB, dx)(2) >= 0.0 && zisa::cross(dx, dC)(2) >= 0.0;
}

TetrahedralCone::TetrahedralCone(const XYZ &A,
                                 const XYZ &B,
                                 const XYZ &C,
                                 const XYZ &D)
    : A(A), dB(XYZ(B - A)), dC(XYZ(C - A)), dD(XYZ(D - A)) {}

bool TetrahedralCone::is_inside(const XYZ &x) const {
  auto dx = XYZ(x - A);

  // clang-format off
  return  zisa::det(dB, dC, dx) >= 0.0
       && zisa::det(dC, dD, dx) >= 0.0
       && zisa::det(dD, dB, dx) >= 0.0;
  // clang-format on
}

bool FullSphere::is_inside(const XYZ &) const { return true; }

Ball::Ball(const XYZ &x_center, double radius)
    : x_center(x_center), radius(radius) {}

bool Ball::is_inside(const XYZ &x) const {
  return zisa::norm(x - x_center) < radius;
}

bool UnionOfRegions::is_inside(const XYZ &x) const {
  for (const auto &r : regions) {
    if (r->is_inside(x)) {
      return true;
    }
  }

  return false;
}

void UnionOfRegions::add(std::shared_ptr<Region> region) {
  regions.push_back(std::move(region));
}
} // namespace zisa
