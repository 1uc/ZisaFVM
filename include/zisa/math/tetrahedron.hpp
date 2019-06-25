#ifndef ZISA_TETRAHEDRON_HPP
#define ZISA_TETRAHEDRON_HPP

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {
class Tetrahedron {
public:
  XYZ points[4];
};

double volume(const Tetrahedron &tetrahedron);
XYZ barycenter(const Tetrahedron &tetrahedron);

}

#endif // ZISA_TETRAHEDRON_HPP
