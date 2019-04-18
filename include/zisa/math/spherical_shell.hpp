#ifndef ZISA_SPHERICAL_SHELL_DDKOWE_HPP
#define ZISA_SPHERICAL_SHELL_DDKOWE_HPP

#include <zisa/config.hpp>

#include <zisa/math/basic_functions.hpp>
#include <zisa/math/mathematical_constants.hpp>

namespace zisa {
struct SphericalShell {
  double r_lower;
  double r_upper;
};

inline double volume(const SphericalShell &shell) {
  auto [r_lower, r_upper] = shell;
  return 4.0 / 3.0 * zisa::pi * (zisa::pow<3>(r_upper) - zisa::pow<3>(r_lower));
}

}

#endif // ZISA_SPHERICAL_SHELL_HPP
