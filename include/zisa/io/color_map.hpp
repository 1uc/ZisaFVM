#ifndef COLOR_MAP_H_XAXAU
#define COLOR_MAP_H_XAXAU

#include <cassert>

#include <zisa/config.hpp>
#include <zisa/io/colors.hpp>

namespace zisa {

class ColorMap {
public:
  ColorMap(int_t n_bins_)
      : n_bins(n_bins_ + (n_bins_ + 1) % 2), colors(n_bins) {

    auto min = rgb2lab(RGBColor{0.0f, 0.353f, 0.498f});
    auto max = rgb2lab(RGBColor{0.498f, 0.165f, 0.192f});

    auto current = min;
    for (int_t i = 0; i < n_bins / 2; ++i) {
      float dL = 100.0f - min.L;
      current.L = min.L + 2.0f * float(i) / float(n_bins) * dL;
      colors[i] = lab2rgb(current);
    }

    colors[n_bins / 2] = lab2rgb(LABColor{100.0f, 0.0f, 0.0f});

    current = max;
    for (int_t i = n_bins / 2 + 1; i < n_bins; ++i) {
      float dL = max.L - 100.f;
      current.L = 100.0f + 2.0f * float(i - n_bins / 2) / float(n_bins) * dL;
      colors[i] = lab2rgb(current);
    }
  }

  /// Color corresponding to the scalar `y`.
  /**
   *  Input:
   *      y : normalized to [-1, 1]
   */
  const RGBColor &operator()(double y) const { return colors[index(y)]; }

private:
  int_t index(double y) const {
    assert(-1.0 <= y);
    assert(y <= 1.0);

    return zisa::min(int_t(0.5 * (y + 1.0) * float(n_bins)), n_bins - 1);
  }

private:
  int_t n_bins;
  std::vector<RGBColor> colors;
};

} // namespace zisa
#endif
