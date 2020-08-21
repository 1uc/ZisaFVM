#ifndef LOCAL_RECONSTRUCTION_H_VF8YB
#define LOCAL_RECONSTRUCTION_H_VF8YB

#include <zisa/grid/grid.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/model/characteristic_scale.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/local_equilibrium.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

// TODO remove
#include <zisa/mpi/mpi.hpp>

namespace zisa {

struct LocalRCParams {
  int_t steps_per_recompute;
  double recompute_threshold;
};

template <class Equilibrium, class RC, class Scaling>
class LocalReconstruction {
private:
  using cvars_t = euler_var_t;

public:
  LocalReconstruction() = default;
  LocalReconstruction(std::shared_ptr<Grid> grid,
                      const LocalEquilibrium<Equilibrium> &eq,
                      const RC &rc,
                      int_t i_cell,
                      Scaling scaling,
                      const LocalRCParams &params)
      : grid(grid),
        eq(eq),
        rc(rc),
        i_cell(i_cell),
        scaling(scaling),
        rhoE_cache(shape_t<1>{rc.local2global().size()}),
        steps_per_recompute(params.steps_per_recompute),
        recompute_threshold(params.recompute_threshold) {}

  void compute_equilibrium(const array_view<cvars_t, 1> &u_local) {
    const auto &u0 = u_local(int_t(0));
    auto rhoE_self = RhoE{u0[0], internal_energy(u0)};

    auto &l2g = rc.local2global();
    scale = scaling(rhoE_self);
    try {
      eq.solve(rhoE_self, grid->cells(i_cell));
    } catch (const std::runtime_error &e) {
      auto msg = string_format("%d %d", zisa::mpi::rank(), l2g[0]);
      PRINT(msg);
      throw e;
    }
    for (int_t il = 0; il < l2g.size(); ++il) {
      try {
        rhoE_cache(il) = eq.extrapolate(grid->cells(l2g[il]));
      } catch (const std::runtime_error &e) {
        auto x0 = grid->cell_centers[l2g[0]];
        auto xil = grid->cell_centers[l2g[il]];
        auto msg1 = string_format("%d %d %s %s",
                                  l2g[0],
                                  l2g[il],
                                  format_as_list(x0).c_str(),
                                  format_as_list(xil).c_str());

        auto rhoE_self_str = format_as_list(rhoE_self);
        auto rhoE_0_str = format_as_list(rhoE_cache(0));
        auto drhoE_str
            = format_as_list(Cartesian<2>(rhoE_self - rhoE_cache(0)));

        auto msg2 = string_format("%s - %s = %s",
                                  rhoE_self_str.c_str(),
                                  rhoE_0_str.c_str(),
                                  drhoE_str.c_str());
        PRINT(msg1);
        PRINT(msg2);
        throw e;
      }
    }

    steps_since_recompute = 0;
  }

  void recompute_equilibrium(const array_view<cvars_t, 1> &u_local) {
    if (steps_since_recompute % steps_per_recompute == 0) {
      compute_equilibrium(u_local);
    } else {
      const auto &u0 = u_local(int_t(0));
      auto [rho, E] = RhoE{u0[0], internal_energy(u0)};
      auto [rho_c, E_c] = rhoE_cache(0);
      auto delta_rhoE = RhoE{(rho - rho_c) / scale[0], (E - E_c) / scale[4]};

      if (zisa::norm(delta_rhoE) >= recompute_threshold) {
        compute_equilibrium(u_local);
      }
    }
  }

  void compute(const array_view<double, 2, row_major> &rhs,
               const array_view<WENOPoly, 1> &polys,
               const array_view<cvars_t, 1> &u_local) {

    recompute_equilibrium(u_local);

    auto &l2g = rc.local2global();
    for (int_t il = 0; il < l2g.size(); ++il) {
      auto [rho_eq_bar, E_eq_bar] = equilibrium_cell_average(il);

      u_local(il)[0] -= rho_eq_bar;
      u_local(il)[4] -= E_eq_bar;

      u_local(il) = u_local(il) / scale;
    }

    weno_poly = rc.reconstruct(rhs, polys, u_local);
    ++steps_since_recompute;
  }

  RhoE equilibrium_cell_average(int_t il) { return rhoE_cache(il); }

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
  RC rc;
  int_t i_cell;
  WENOPoly weno_poly;
  array<ScalarPoly, 1> scalar_polys;
  Scaling scaling;
  cvars_t scale = cvars_t::zeros();

  array<RhoE, 1> rhoE_cache;
  int_t steps_per_recompute;
  int_t steps_since_recompute = 0;
  double recompute_threshold;
};

} // namespace zisa

#endif /* end of include guard */
