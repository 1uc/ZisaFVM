#ifndef ZISA_TETRAHEDRON_HPP
#define ZISA_TETRAHEDRON_HPP

#include <zisa/config.hpp>
#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/triangle.hpp>

namespace zisa {
class Tetrahedron {
public:
  XYZ points[4];

  Tetrahedron(const XYZ &v0, const XYZ &v1, const XYZ &v2, const XYZ &v3);
};

double inradius(const Tetrahedron &tet);

Triangle face(const Tetrahedron &tet, int_t k);
double volume(const Tetrahedron &tetrahedron);

XYZ barycenter(const Tetrahedron &tetrahedron);

double characteristic_length(const Tetrahedron &tet);

}

#endif // ZISA_TETRAHEDRON_HPP
