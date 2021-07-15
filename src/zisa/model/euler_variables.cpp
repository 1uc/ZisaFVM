// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#include <zisa/math/face.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

std::string euler_var_t::labels(int_t i) {
  if (i == 0) {
    return "rho";
  }
  if (i == 1) {
    return "mv1";
  }
  if (i == 2) {
    return "mv2";
  }
  if (i == 3) {
    return "mv3";
  }
  if (i == 4) {
    return "E";
  }

  LOG_ERR(string_format("Unknown input. [%d]", i));
}

void coord_transform(euler_var_t &u,
                     const XYZ &n,
                     const XYZ &t1,
                     const XYZ &t2) {

  double u_dot_n = u[1] * n[0] + u[2] * n[1] + u[3] * n[2];
  double u_dot_t1 = u[1] * t1[0] + u[2] * t1[1] + u[3] * t1[2];
  double u_dot_t2 = u[1] * t2[0] + u[2] * t2[1] + u[3] * t2[2];

  u[1] = u_dot_n;
  u[2] = u_dot_t1;
  u[3] = u_dot_t2;
}

void coord_transform(euler_var_t &u, const Face &face) {
  const auto &n = face.normal;
  const auto &[t1, t2] = face.tangentials;
  coord_transform(u, n, t1, t2);
}

void inv_coord_transform(euler_var_t &u,
                         const XYZ &n,
                         const XYZ &t1,
                         const XYZ &t2) {
  double ux = u[1] * n[0] + u[2] * t1[0] + u[3] * t2[0];
  double uy = u[1] * n[1] + u[2] * t1[1] + u[3] * t2[1];
  double uz = u[1] * n[2] + u[2] * t1[2] + u[3] * t2[2];

  u[1] = ux;
  u[2] = uy;
  u[3] = uz;
}

void inv_coord_transform(euler_var_t &u, const Face &face) {
  const auto &n = face.normal;
  const auto &[t1, t2] = face.tangentials;

  inv_coord_transform(u, n, t1, t2);
}

} // namespace zisa
