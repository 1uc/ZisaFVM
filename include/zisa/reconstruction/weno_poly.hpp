#ifndef WENO_POLY_H_WYSZQ
#define WENO_POLY_H_WYSZQ

#include <zisa/config.hpp>

#include <zisa/math/poly2d.hpp>

namespace zisa {
using WENOPoly = PolyND</* max_coeffs = */ 35, /* vars = */ 5>;
using ScalarPoly = PolyND</* max_coeffs = */ 35, /* vars = */ 1>;
}

#endif
