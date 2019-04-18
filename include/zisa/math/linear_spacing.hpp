#ifndef ZISA_LINEAR_SPACING_HPP
#define ZISA_LINEAR_SPACING_HPP

#include <zisa/config.hpp>
#include <zisa/memory/array.hpp>

namespace zisa {

array<double, 1> linear_spacing(double x0, double x1, int_t n_points);

}

#endif // ZISA_LINEAR_SPACING_HPP
