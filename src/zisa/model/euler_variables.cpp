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

void coord_transform(euler_var_t &u, const XYZ &n, const XYZ &t) {
  double u_dot_n = u[1] * n[0] + u[2] * n[1];
  double u_dot_t = u[1] * t[0] + u[2] * t[1];

  u[1] = u_dot_n;
  u[2] = u_dot_t;
}

void coord_transform(euler_var_t &u, const Edge &edge) {
  coord_transform(u, edge.normal(), edge.tangential());
}

void inv_coord_transform(euler_var_t &u, const XYZ &n, const XYZ &t) {
  double ux = u[1] * n[0] + u[2] * t[0];
  double uy = u[1] * n[1] + u[2] * t[1];

  u[1] = ux;
  u[2] = uy;
}

void inv_coord_transform(euler_var_t &u, const Edge &edge) {
  inv_coord_transform(u, edge.normal(), edge.tangential());
}

} // namespace zisa
