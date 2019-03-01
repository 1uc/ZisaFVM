#ifndef ZISA_BOUNDING_BOX_HPP
#define ZISA_BOUNDING_BOX_HPP

#include <zisa/grid/grid.hpp>
#include <zisa/math/cartesian.hpp>

namespace zisa {

struct BoundingBox {
  XYZ min;
  XYZ max;
};

BoundingBox bounding_box(const Grid &grid);

}
#endif // ZISA_BOUNDING_BOX_HPP
