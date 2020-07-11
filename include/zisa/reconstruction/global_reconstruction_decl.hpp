#ifndef GLOBAL_RECONSTRUCTION_DECL_H_VHVV5
#define GLOBAL_RECONSTRUCTION_DECL_H_VHVV5

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/memory/block_allocator.hpp>
#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/local_eos_state.hpp>
#include <zisa/reconstruction/hybrid_weno_params.hpp>
#include <zisa/reconstruction/local_reconstruction.hpp>
#include <zisa/reconstruction/stencil_family.hpp>
#include <zisa/reconstruction/weno_poly.hpp>

namespace zisa {

template <class CVars>
class GlobalReconstruction {
public:
  virtual ~GlobalReconstruction() = default;

  virtual void compute(const AllVariables &current_state) = 0;
  virtual CVars operator()(int_t i, const XYZ &x) const = 0;
};

template <class Equilibrium, class RC, class Scaling>
class EulerGlobalReconstruction : public GlobalReconstruction<euler_var_t> {
private:
  using cvars_t = euler_var_t;
  using lrc_t = LocalReconstruction<Equilibrium, RC, Scaling>;

public:
  EulerGlobalReconstruction(const HybridWENOParams &params, array<lrc_t, 1> rc);

  const lrc_t &operator()(int_t i) const;

  virtual euler_var_t operator()(int_t i, const XYZ &x) const override;
  virtual void compute(const AllVariables &current_state) override;

  std::string str() const;

private:
  void set_qbar_local(array<cvars_t, 1> &qbar_local,
                      const AllVariables &current_state,
                      int_t i);

  void set_tracer_local(array<double, 2, column_major> &tracer_local,
                        const AllVariables &current_state,
                        int_t i);

private:
  HybridWENOParams params;
  int_t max_stencil_size;
  int_t n_polys;

  array<lrc_t, 1> rc;
  std::shared_ptr<block_allocator<array<cvars_t, 1>>> qbar_allocator;
  std::shared_ptr<block_allocator<array<double, 2, column_major>>>
      tracer_allocator;
  std::shared_ptr<block_allocator<array<WENOPoly, 1>>> polys_allocator;
  std::shared_ptr<block_allocator<array<double, 2, row_major>>> rhs_allocator;
};

template <class Equilibrium, class RC, class Scaling, class EOS, class Gravity>
array<LocalReconstruction<Equilibrium, RC, Scaling>, 1>
make_reconstruction_array(const std::shared_ptr<Grid> &grid,
                          const HybridWENOParams &weno_params,
                          const LocalEOSState<EOS> &local_eos_state,
                          const std::shared_ptr<Gravity> &gravity) {

  auto n_cells = grid->n_cells;

  auto lrc = array<LocalReconstruction<Equilibrium, RC, Scaling>, 1>(
      shape_t<1>{n_cells});
  for (int_t i_cell = 0; i_cell < n_cells; ++i_cell) {
    const auto &eos = local_eos_state(i_cell);

    auto rc = RC(grid, i_cell, weno_params);
    lrc[i_cell] = LocalReconstruction<Equilibrium, RC, Scaling>(
        grid,
        LocalEquilibrium<Equilibrium>(
            make_equilibrium<Equilibrium>(eos, gravity)),
        rc,
        i_cell,
        Scaling(eos));
  }

  return lrc;
}

template <class Equilibrium, class RC, class Scaling, class EOS, class Gravity>
array<LocalReconstruction<Equilibrium, RC, Scaling>, 1>
make_reconstruction_array(const std::shared_ptr<Grid> &grid,
                          const array<StencilFamily, 1> &stencils,
                          const HybridWENOParams &weno_params,
                          const LocalEOSState<EOS> &local_eos_state,
                          const std::shared_ptr<Gravity> &gravity) {

  auto n_cells = grid->n_cells;

  auto lrc = array<LocalReconstruction<Equilibrium, RC, Scaling>, 1>(
      shape_t<1>{n_cells});
  for (int_t i = 0; i < n_cells; ++i) {
    const auto &eos = local_eos_state(i);
    auto eq = make_equilibrium<Equilibrium>(eos, gravity);
    auto scaling = Scaling(eos);

    auto o1_params = HybridWENOParams(
        {{1}, {"c"}, {1.0}}, {1.0}, weno_params.epsilon, weno_params.exponent);

    if (stencils[i].size() == 1 && stencils[i].order() == 1) {
      lrc[i] = LocalReconstruction<Equilibrium, RC, Scaling>(
          grid,
          LocalEquilibrium(eq),
          RC(grid, stencils[i], i, o1_params),
          i,
          scaling);

    } else {
      lrc[i] = LocalReconstruction<Equilibrium, RC, Scaling>(
          grid,
          LocalEquilibrium(eq),
          RC(grid, stencils[i], i, weno_params),
          i,
          scaling);
    }
  }

  return lrc;
}

} // namespace zisa
#endif
