#ifndef COLORS_H_8XS7V
#define COLORS_H_8XS7V

#include <array>
#include <ostream>

#include <zisa/config.hpp>
#include <zisa/math/basic_functions.hpp>

namespace zisa {

struct SRGBColor {
  float r;
  float g;
  float b;
};

std::ostream &operator<<(std::ostream &os, const SRGBColor &srgb);

struct RGBColor {
  float r;
  float g;
  float b;
};

std::ostream &operator<<(std::ostream &os, const RGBColor &rgb);
bool is_valid(const RGBColor &rgb);
bool almost_equal(const RGBColor &lhs, const RGBColor &rhs, float atol);

struct XYZColor {
  float x;
  float y;
  float z;
};

std::ostream &operator<<(std::ostream &os, const XYZColor &xyz);
bool is_valid(const XYZColor &xyz);
bool almost_equal(const XYZColor &lhs, const XYZColor &rhs, float atol);

struct LABColor {
  float L;
  float a;
  float b;
};

std::ostream &operator<<(std::ostream &os, const LABColor &lab);
bool is_valid(const LABColor &lab);
bool almost_equal(const LABColor &lhs, const LABColor &rhs, float atol);

// RGBColor srgb2rgb(const SRGBColor &srgb);
// SRGBColor rgb2srgb(const RGBColor &rgb);

XYZColor rgb2xyz(const RGBColor &rgb);
RGBColor xyz2rgb(const XYZColor &xyz);

LABColor xyz2lab(const XYZColor &xyz);
XYZColor lab2xyz(const LABColor &lab);

LABColor rgb2lab(const RGBColor &rgb);
RGBColor lab2rgb(const LABColor &lab);

} // namespace zisa

#endif /* end of include guard */
