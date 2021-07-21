// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/io/colors.hpp>

namespace zisa {

std::ostream &operator<<(std::ostream &os, const SRGBColor &srgb) {
  os << string_format("[%e, %e, %e]", srgb.r, srgb.g, srgb.b);
  return os;
}

std::ostream &operator<<(std::ostream &os, const RGBColor &rgb) {
  os << string_format("[%e, %e, %e]", rgb.r, rgb.g, rgb.b);
  return os;
}

bool is_valid(const RGBColor &rgb) {
  float r = rgb.r, g = rgb.g, b = rgb.b;

  return (0.0f <= r && r <= 1.0f) && (0.0f <= g && g <= 1.0f)
         && (0.0f <= b && b <= 1.0f);
}

bool almost_equal(const RGBColor &lhs, const RGBColor &rhs, float atol) {
  return almost_equal(lhs.r, rhs.r, atol) && almost_equal(lhs.g, rhs.g, atol)
         && almost_equal(lhs.b, rhs.b, atol);
}

std::ostream &operator<<(std::ostream &os, const XYZColor &xyz) {
  os << string_format("[%e, %e, %e]", xyz.z, xyz.y, xyz.z);
  return os;
}

bool almost_equal(const XYZColor &lhs, const XYZColor &rhs, float atol) {
  return almost_equal(lhs.x, rhs.x, atol) && almost_equal(lhs.y, rhs.y, atol)
         && almost_equal(lhs.z, rhs.z, atol);
}

bool is_valid(const XYZColor &xyz) {
  float x = xyz.x, y = xyz.y, z = xyz.z;

  return 0.0f < x && 0.0f < y && 0.0f < z;
}

std::ostream &operator<<(std::ostream &os, const LABColor &lab) {
  os << string_format("[%e, %e, %e]", lab.L, lab.a, lab.b);
  return os;
}

bool is_valid(const LABColor &lab) {
  float L = lab.L, a = lab.a, b = lab.b;

  return (0.0f <= L && L <= 100.0f) && (-100.0f <= a && a <= 100.0f)
         && (-100.0f <= b && b <= 100.0f);
}

bool almost_equal(const LABColor &lhs, const LABColor &rhs, float atol) {
  return almost_equal(lhs.L, rhs.L, atol) && almost_equal(lhs.a, rhs.a, atol)
         && almost_equal(lhs.b, rhs.b, atol);
}

// RGBColor srgb2rgb(const SRGBColor &srgb);
// SRGBColor rgb2srgb(const RGBColor &rgb);

XYZColor rgb2xyz(const RGBColor &rgb) {
  /// Reference: Kenneth Moreland, Sandia National Labs

  float r = rgb.r, g = rgb.g, b = rgb.b;

  float x = 100.0f * (0.4124f * r + 0.3576f * g + 0.1805f * b);
  float y = 100.0f * (0.2126f * r + 0.7152f * g + 0.0722f * b);
  float z = 100.0f * (0.0193f * r + 0.1192f * g + 0.9505f * b);

  return XYZColor{x, y, z};
}

RGBColor xyz2rgb(const XYZColor &xyz) {
  float x = xyz.x, y = xyz.y, z = xyz.z;

  float r = 0.01f * (3.24062548f * x + -1.53720797f * y + -0.4986286f * z);
  float g = 0.01f * (-0.96893071f * x + 1.87575606f * y + 0.04151752f * z);
  float b = 0.01f * (0.05571012f * x + -0.20402105f * y + 1.05699594f * z);

  return RGBColor{r, g, b};
}

namespace cielab {

std::array<float, 3> reference_white() { return {95.047f, 100.0f, 108.883f}; }

float f(float chi) {
  return chi > 0.008856f ? zisa::pow(chi, 1.0f / 3.0f)
                         : 7.787f * chi + 16.0f / 116.0f;
}

float finv(float chi) {
  return chi < 0.206892f ? (chi - 16.0f / 116.0f) / 7.787f : chi * chi * chi;
}

} // namespace cielab

LABColor xyz2lab(const XYZColor &xyz) {
  const auto &[xn, yn, zn] = cielab::reference_white();

  float x = xyz.x, y = xyz.y, z = xyz.z;

  auto fy = cielab::f(y / yn);

  float L = 116.0f * (fy - 16.0f / 116.0f);
  float a = 500.0f * (cielab::f(x / xn) - fy);
  float b = 200.0f * (fy - cielab::f(z / zn));

  return LABColor{L, a, b};
}

XYZColor lab2xyz(const LABColor &lab) {
  const auto &[xn, yn, zn] = cielab::reference_white();

  float L = lab.L, a = lab.a, b = lab.b;

  float y = yn * cielab::finv(L / 116.0f + 16.0f / 116.0f);
  float fy = cielab::f(y / yn);

  float x = xn * cielab::finv(a / 500.0f + fy);
  float z = zn * cielab::finv(-b / 200.f + fy);

  return XYZColor{x, y, z};
}

LABColor rgb2lab(const RGBColor &rgb) { return xyz2lab(rgb2xyz(rgb)); }
RGBColor lab2rgb(const LABColor &lab) { return xyz2rgb(lab2xyz(lab)); }

} // namespace zisa
