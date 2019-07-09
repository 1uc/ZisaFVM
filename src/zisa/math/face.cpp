#include <zisa/math/face.hpp>

#include <zisa/math/quadrature.hpp>

namespace zisa {
XYZ barycenter(const Face &face) {
  return average(face.qr, [](const XYZ &x) { return x; });
}

Face::Face(DenormalizedRule qr, const XYZ &n, const XYZ &t1, const XYZ &t2)
    : qr(std::move(qr)), normal(n), tangentials{t1, t2} {}

bool operator==(const Face &a, const Face &b) {
  return a.qr == b.qr && a.normal == b.normal && a.tangentials == b.tangentials;
}

bool operator!=(const Face &a, const Face &b) { return !(a == b); }

double volume(const Face &face) { return volume(face.qr); }

XYZ unit_outward_normal(const Face &face, const XYZ &point_inside) {
  const auto &n = face.normal;
  return XYZ(zisa::sign(zisa::dot(n, face.qr.points[0] - point_inside)) * n);
}
}