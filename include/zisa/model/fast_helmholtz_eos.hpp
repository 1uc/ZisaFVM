// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_FAST_HELMHOLTZ_EOS_HPP_CCCU
#define ZISA_FAST_HELMHOLTZ_EOS_HPP_CCCU

#if ZISA_HAS_HELMHOLTZ_EOS != 1
#error "Missing `ZISA_HAS_HELMHOLTZ_EOS=1`."
#endif

#include <zisa/config.hpp>
#include <zisa/model/helmholtz_eos.hpp>

namespace zisa {

class FastHelmholtzEOS : public EquationOfState {
private:
  using super = EquationOfState;

public:
  using cvars_t = super::cvars_t;
  using xvars_t = super::xvars_t;



private:
  double a_bar;
  double z_bar;

  double gamma_p;
  double gamma_cs;
};

}
#endif // ZISA_FAST_HELMHOLTZ_EOS_HPP