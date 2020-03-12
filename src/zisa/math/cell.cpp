#include <zisa/math/cell.hpp>

#include <zisa/math/quadrature.hpp>

namespace zisa {
XYZ barycenter(const Cell &cell) {
  return average(cell, [](const XYZ &x) { return x; });
}

double avg_moment(const Cell &cell, int x_deg, int y_deg) {
  return avg_moment(cell, x_deg, y_deg, 0);
}

double avg_moment(const Cell &cell, int x_deg, int y_deg, int z_deg) {
  auto center = barycenter(cell);

  auto f = [x_deg, y_deg, z_deg, &center](const XYZ &xyz) {
    auto [x, y, z] = XYZ(xyz - center);
    return zisa::pow(x, x_deg) * zisa::pow(y, y_deg) * zisa::pow(z, z_deg);
  };

  return average(cell, f);
}

std::string str(const Cell &cell) {
    return string_format("qr = {%s}", str(cell.qr).c_str());
}

}
