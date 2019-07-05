#include <zisa/math/tetrahedron.hpp>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/loops/execution_policies.hpp>
#include <zisa/loops/reduction/min.hpp>

namespace zisa {
Tetrahedron::Tetrahedron(const XYZ &v0,
                         const XYZ &v1,
                         const XYZ &v2,
                         const XYZ &v3)
    : points{v0, v1, v2, v3} {}

double volume(const Tetrahedron &tetrahedron) {
  const auto &[v0, v1, v2, v3] = tetrahedron.points;

  auto [d11, d12, d13] = XYZ(v1 - v0);
  auto [d21, d22, d23] = XYZ(v2 - v0);
  auto [d31, d32, d33] = XYZ(v3 - v0);

  // clang-format off
  double vol = 1.0/6.0 * (d11*d22*d33 + d21*d32*d13 + d31*d12*d23
                        - d11*d32*d23 - d21*d12*d33 - d31*d22*d13);
  // clang-format on

  return zisa::abs(vol);
}

XYZ barycenter(const Tetrahedron &tetrahedron) {
  const auto &[v0, v1, v2, v3] = tetrahedron.points;
  return XYZ(0.25 * (v0 + v1 + v2 + v3));
}

double characteristic_length(const Tetrahedron &tet) {
  auto c = barycenter(tet);

  auto r = [&c](const XYZ &x) { return 2.0 * zisa::norm(c - x); };

  double l = 0.0;
  for (const auto &p : tet.points) {
    l += r(p);
  }

  return 0.25 * l;
}

double inradius(const Tetrahedron &tet) {
  auto c = barycenter(tet);

  return zisa::reduce::min(
      serial_policy{}, PlainIndexRange(0, 4), [&c, &tet](int_t k) {
        return zisa::norm(barycenter(face(tet, k)) - c);
      });
}

Triangle face(const Tetrahedron &tet, int_t k) {
  auto element_type = GMSHElementType::tetrahedron;
  auto iv0 = GMSHElementInfo::relative_vertex_index(element_type, k, 0);
  auto iv1 = GMSHElementInfo::relative_vertex_index(element_type, k, 1);
  auto iv2 = GMSHElementInfo::relative_vertex_index(element_type, k, 2);

  const auto &v0 = tet.points[iv0];
  const auto &v1 = tet.points[iv1];
  const auto &v2 = tet.points[iv2];

  return Triangle(v0, v1, v2);
}

}
