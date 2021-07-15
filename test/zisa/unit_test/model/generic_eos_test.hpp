// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_GENERIC_EOS_TEST_XIMKL_HPP
#define ZISA_GENERIC_EOS_TEST_XIMKL_HPP

#include <zisa/model/euler_variables.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/utils/to_string.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

namespace detail {
template <class To>
class Convert;

template <>
struct Convert<RhoP> {
  template <class EOS, class From>
  static RhoP convert(const EOS &eos, const From &from) {
    return eos.rhoP(from);
  }
};

template <>
struct Convert<RhoEntropy> {
  template <class EOS, class From>
  static RhoEntropy convert(const EOS &eos, const From &from) {
    return eos.rhoK(from);
  }
};

template <>
struct Convert<RhoE> {
  template <class EOS, class From>
  static RhoE convert(const EOS &eos, const From &from) {
    return eos.rhoE(from);
  }
};

template <>
struct Convert<EnthalpyEntropy> {
  template <class EOS, class From>
  static EnthalpyEntropy convert(const EOS &eos, const From &from) {
    return eos.enthalpy_entropy(from);
  }
};
}

template <class To, class EOS, class From>
To convert(const EOS &eos, const From &from) {
  return detail::Convert<To>::convert(eos, from);
}

}

template <class To, class EOS, class From>
void check_eos(const EOS &eos, const From &exact, double rtol = 1e-10) {
  To to = zisa::convert<To>(eos, exact);
  From approx = zisa::convert<From>(eos, to);

  auto decl = string_format("%s -> %s -> %s",
                            type_name<From>().c_str(),
                            type_name<To>().c_str(),
                            type_name<From>().c_str());

  auto values = string_format("%s -> %s -> %s  err = %s",
                              zisa::to_string(exact).c_str(),
                              zisa::to_string(to).c_str(),
                              zisa::to_string(approx).c_str(),
                              zisa::to_string(From(exact - approx)).c_str());

  INFO(string_format(
      "%s\n%s\n%s", type_name<EOS>().c_str(), decl.c_str(), values.c_str()));
  REQUIRE(zisa::almost_equal(approx, exact, rtol * zisa::norm(exact)));
}

template <class EOS>
void generic_eos_tests(
    const EOS &eos, double rho, double p, double E, double K, double h) {

  double rtol = 1e-8;

  check_eos<zisa::RhoP>(eos, zisa::RhoE{rho, E});
  check_eos<zisa::RhoP>(eos, zisa::RhoEntropy{rho, K});
  check_eos<zisa::RhoP>(eos, zisa::EnthalpyEntropy{h, K}, rtol);

  check_eos<zisa::RhoE>(eos, zisa::RhoP{rho, p});
  check_eos<zisa::RhoE>(eos, zisa::RhoEntropy{rho, K});
  check_eos<zisa::RhoE>(eos, zisa::EnthalpyEntropy{h, K}, rtol);

  check_eos<zisa::RhoEntropy>(eos, zisa::RhoP{rho, p});
  check_eos<zisa::RhoEntropy>(eos, zisa::RhoE{rho, E});
  check_eos<zisa::RhoEntropy>(eos, zisa::EnthalpyEntropy{h, K}, rtol);

  check_eos<zisa::EnthalpyEntropy>(eos, zisa::RhoP{rho, p}, rtol);
  check_eos<zisa::EnthalpyEntropy>(eos, zisa::RhoE{rho, E}, rtol);
  check_eos<zisa::EnthalpyEntropy>(eos, zisa::RhoEntropy{rho, K}, rtol);
}

#endif // ZISA_GENERIC_EOS_TEST_HPP
