#ifndef POLYTROPE_IC_H_PK11F
#define POLYTROPE_IC_H_PK11F

#include <zisa/math/cartesian.hpp>
#include <zisa/model/euler_variables.hpp>

namespace zisa {

template <class EULER>
class PolytropeIC {
private:
  using euler_t = EULER;

public:
  PolytropeIC(const euler_t &euler) : euler(euler) {}

  RhoP operator()(const XYZ &x) const {

    double alpha = this->euler.gravity.alpha();
    double eps = std::numeric_limits<double>::min();
    double r_eff = alpha * (zisa::norm(x) + eps);

    double rho = zisa::sin(r_eff) / r_eff;
    double p = zisa::pow<2>(rho);

    return RhoP{rho, p};
  }

private:
  euler_t euler;
};

} // namespace zisa
#endif /* end of include guard */
