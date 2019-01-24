#ifndef WENO_POLY_H_WYSZQ
#define WENO_POLY_H_WYSZQ

#include <zisa/config.hpp>

#include <zisa/math/poly2d.hpp>

namespace zisa {
using WENOPoly = Poly2D</* degree = */ 4, /* vars = */ 5>;
}

#endif
