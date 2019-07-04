#ifndef LOCAL_RECONSTRUCTION_H_VF8YB
#define LOCAL_RECONSTRUCTION_H_VF8YB

#include <zisa/grid/grid.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/characteristic_scale.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/local_equilibrium.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {

template <class Equilibrium, class RC, class Scaling>
class LocalReconstruction {
private:
  using cvars_t = euler_var_t;

public:
  LocalReconstruction() = default;
  LocalReconstruction(std::shared_ptr<Grid> &grid,
                      const LocalEquilibrium<Equilibrium> &eq,
                      const RC &rc,
                      const Triangle &tri_ref,
                      Scaling scaling)
      : grid(grid), eq(eq), rc(rc), tri_ref(tri_ref), scaling(scaling) {}

  void compute(array<cvars_t, 1> &u_local) {
    const auto &u0 = u_local(int_t(0));
    auto rhoE_self = RhoE{u0[0], internal_energy(u0)};
    eq.solve(rhoE_self, tri_ref);

    scale = scaling(rhoE_self);

    auto &l2g = rc.local2global();
    for (int_t il = 0; il < l2g.size(); ++il) {
      auto [rho_eq_bar, E_eq_bar] = eq.extrapolate(triangle(*grid, l2g[il]));

      u_local(il)[0] -= rho_eq_bar;
      u_local(il)[4] -= E_eq_bar;

      u_local(il) = u_local(il) / scale;
    }

    weno_poly = rc.reconstruct(u_local);
  }

  cvars_t operator()(const XYZ &x) const {
    return cvars_t(background(x) + delta(x));
  }

  cvars_t delta(const XYZ &x) const { return cvars_t(scale * weno_poly(x)); }

  cvars_t background(const XYZ &x) const {
    auto [rho, E] = eq.extrapolate(x);
    return cvars_t{rho, 0.0, 0.0, 0.0, E};
  }

  auto combined_stencil_size() const
      -> decltype(std::declval<RC>().combined_stencil_size()) {
    return rc.combined_stencil_size();
  }

  auto local2global() const -> decltype(std::declval<RC>().local2global()) {
    return rc.local2global();
  }

private:
  std::shared_ptr<Grid> grid;
  LocalEquilibrium<Equilibrium> eq;
  RC rc;
  Triangle tri_ref;
  WENOPoly weno_poly;
  Scaling scaling;
  cvars_t scale = cvars_t::zeros();
};

} // namespace zisa

#endif /* end of include guard */
