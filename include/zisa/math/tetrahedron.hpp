#ifndef ZISA_TETRAHEDRON_HPP
#define ZISA_TETRAHEDRON_HPP

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {
class Tetrahedron {
public:
  XYZ points[4];

  Tetrahedron(const XYZ &v0, const XYZ &v1, const XYZ &v2, const XYZ &v3);
};

double volume(const Tetrahedron &tetrahedron);
XYZ barycenter(const Tetrahedron &tetrahedron);

double characteristic_length(const Tetrahedron &tet);

}

#endif // ZISA_TETRAHEDRON_HPP
