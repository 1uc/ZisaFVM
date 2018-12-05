#include <zisa/io/scalar_plot.hpp>

namespace zisa {

ScalarPlot::ScalarPlot(const array<XY, 1> &vertices,
                       const array<int_t, 2> &vertex_indices,
                       const std::string &title)
    : tri_plot(vertices, vertex_indices, title),
      color_map(256),
      values(shape_t<1>{vertex_indices.shape(0)}),
      colors(shape_t<1>{3*vertex_indices.shape(0)}) {}

std::pair<double, double> ScalarPlot::safe_min_max(double min_, double max_) {

  // wlog min_ is not denormal
  if (!std::isnormal(min_) && std::isfinite(min_)) {
    min_ = 0.0;
  }

  // wlog max_ is not denormal
  if (!std::isnormal(max_) && std::isfinite(max_)) {
    max_ = 0.0;
  }

  if (min_ != max_) {
    return {min_, max_};
  }

  if (min_ == 0.0) {
    return {-std::numeric_limits<double>::min(),
            std::numeric_limits<double>::min()};
  }

  double eps = std::numeric_limits<double>::epsilon();
  double dx
      = zisa::max(eps * zisa::abs(min_), std::numeric_limits<double>::min());

  return {min_ - dx, max_ + dx};
}

} // namespace zisa
