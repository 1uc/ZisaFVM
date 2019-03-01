#include <zisa/math/bounding_box.hpp>

namespace zisa {

BoundingBox bounding_box(const Grid &grid) {

  auto lower = XYZ{std::numeric_limits<double>::max(),
                   std::numeric_limits<double>::max(),
                   std::numeric_limits<double>::max()};

  auto upper = XYZ{std::numeric_limits<double>::lowest(),
                   std::numeric_limits<double>::lowest(),
                   std::numeric_limits<double>::lowest()};

  for (const auto &v : grid.vertices) {
    for (int_t k = 0; k < XYZ::size(); ++k) {
      lower[k] = zisa::min(v[k], lower[k]);
      upper[k] = zisa::max(v[k], upper[k]);
    }
  }

  return BoundingBox{lower, upper};
}

}