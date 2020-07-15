#ifndef ZISA_CACHED_LOCAL_RECONSTRUCTION_HPP_ZOQENI
#define ZISA_CACHED_LOCAL_RECONSTRUCTION_HPP_ZOQENI

#include <zisa/grid/grid.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/characteristic_scale.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/local_equilibrium.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {

template <class Equilibrium, class RC, class Scaling>
class CachedLocalReconstruction {
private:
  using cvars_t = euler_var_t;

public:
  CachedLocalReconstruction() = default;
  CachedLocalReconstruction(std::shared_ptr<Grid> grid,
                      const LocalEquilibrium<Equilibrium> &eq,
                      const RC &rc,
                      int_t i_cell,
                      Scaling scaling)
      : grid(grid), eq(eq), rc(rc), i_cell(i_cell), scaling(scaling) {

    u_eq_bar = array<cvars_t, 1>(rc.l2g().size());
  }

  void compute(const array_view<double, 2, row_major> &rhs,
               const array_view<WENOPoly, 1> &polys,
               const array_view<cvars_t, 1> &u_local) {

    if (steps_since_recompute % steps_per_recompute == 0) {
            recompute_cache(rhs, polys, u_local);
    } else {
      const auto &u0 = u_local(int_t(0));
      auto [rho, E] = RhoE{u0[0], internal_energy(u0)};
      auto [rho_c, E_c] = rhoE_eq_cache(0);
      auto elta_rhoE = RhoE{(rho - rho_c) / scale[0],
                        (E - E_c) / scale[4]};

      if (zisa::norm(delta_rhoE) >= recompute_threshold) {
        recompute_cache(rhs, polys, u_local);
      }
    }

    auto &l2g = rc.local2global();
    for (int_t il = 0; il < l2g.size(); ++il) {
      auto [rho_eq_bar, E_eq_bar] = rhoE_eq_cache(il);

      u_local(il)[0] -= rho_eq_bar;
      u_local(il)[4] -= E_eq_bar;

      u_local(il) = u_local(il) / scale;
    }

    weno_poly = rc.reconstruct(rhs, polys, u_local);

    ++steps_since_recomputes;
  }


  void recompute_cache(const array_view<double, 2, row_major> &rhs,
               const array_view<WENOPoly, 1> &polys,
               const array_view<cvars_t, 1> &u_local) {

    const auto &u0 = u_local(int_t(0));
    auto rhoE_self = RhoE{u0[0], internal_energy(u0)};

    scale = scaling(rhoE_self);
    eq.solve(rhoE_self, grid->cells(i_cell));
    auto &l2g = rc.local2global();
    for (int_t il = 0; il < l2g.size(); ++il) {
      rhoE_eq_cache(il) = eq.extrapolate(grid->cells(l2g[il]));
    }

    steps_since_recompute = 0;
  }

  void compute_tracer(const array_view<double, 2, row_major> &rhs,
                      const array_view<ScalarPoly, 1> &polys,
                      const array_view<double, 2, column_major> &q_local) {

    auto n_vars = q_local.shape(1);
    if (scalar_polys.size() != n_vars) {
      scalar_polys = array<ScalarPoly, 1>(n_vars);
    }

    auto rhs_view = [&]() {
      auto shape = shape_t<2>{rhs.shape(0), 1};
      return array_view<double, 2, row_major>(shape, rhs.raw());
    }();

    for (int_t k_var = 0; k_var < n_vars; ++k_var) {
      auto q_component = array_const_view<double, 1>(
          shape_t<1>(q_local.shape(0)),
          q_local.raw() + k_var * q_local.shape(0));
      scalar_polys[k_var] = rc.reconstruct(rhs_view, polys, q_component);
    }
  }

  cvars_t operator()(const XYZ &x) const {
    return cvars_t(background(x) + delta(x));
  }

  double tracer(const XYZ &x, int_t k_var) const {
    return scalar_polys[k_var](x)[0];
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

  std::string str(int verbose = 0) const {
    std::stringstream ss;

    ss << grid->str() << "\n";
    ss << eq.str(verbose) << "\n";
    ss << rc.str(verbose) << "\n";

    return ss.str();
  }

private:
  std::shared_ptr<Grid> grid;
  LocalEquilibrium<Equilibrium> eq;

  array<cvars_t, 1> rhoE_eq_cache;
  int_t steps_since_recompute;
  int_t steps_per_recompute;
  double recompute_threshold;

  RC rc;
  int_t i_cell;
  WENOPoly weno_poly;
  array<ScalarPoly, 1> scalar_polys;
  Scaling scaling;
  cvars_t scale = cvars_t::zeros();
};

} // namespace zisa

#endif /* end of include guard */
