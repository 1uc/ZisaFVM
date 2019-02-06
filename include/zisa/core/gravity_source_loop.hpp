#ifndef GRAVITY_SOURCE_LOOP_H_F39RG
#define GRAVITY_SOURCE_LOOP_H_F39RG

#include <zisa/model/all_variables.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {

template <class Equilibrium, class RC, class EULER>
class GravitySourceLoop : public RateOfChange {
private:
  using euler_t = EULER;
  using cvars_t = euler_var_t;
  using grc_t = GlobalReconstruction<Equilibrium, RC>;

public:
  GravitySourceLoop(std::shared_ptr<Grid> grid,
                    const euler_t &euler,
                    std::shared_ptr<grc_t> global_reconstruction,
                    int_t edge_deg,
                    int_t volume_deg)
      : euler(euler),
        grid(grid),
        global_reconstruction(global_reconstruction),
        edge_deg(edge_deg),
        volume_deg(volume_deg) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables & /* current_state */,
                       double /* t */) const override {

    const auto &eos = euler.eos;
    const auto &gravity = euler.gravity;

#pragma omp parallel for schedule(guided)
    for (int_t i = 0; i < grid->n_cells; ++i) {
      const auto &tri = grid->triangle(i);
      const auto &rc = (*global_reconstruction)(i);

      auto x_cell = grid->cell_centers(i);

      // equilibrium terms
      auto s = cvars_t(0.0);
      for (int_t k = 0; k < grid->max_neighbours; ++k) {
        auto edge = grid->edge(i, k);

        auto s_eq = [&x_cell, &rc, &eos, &edge = edge](XY x) {
          auto u_eq = rc.background(x);
          auto p_eq = eos.pressure(RhoE{u_eq[0], u_eq[4]});

          auto n = unit_outward_normal(edge, x_cell);
          auto s = cvars_t{0.0, p_eq * n[0], p_eq * n[1], 0.0, 0.0};

          return s;
        };

        s += quadrature(s_eq, edge, edge_deg);
      }

      // delta terms
      auto s_delta = [&rc, &gravity](XY x) {
        static_assert(XY::size() == 2);

        auto u_eq = rc.background(x);
        auto du = rc.delta(x);
        auto u = cvars_t(u_eq + du);

        auto drho = du[0];
        auto mv = momentum(u);
        auto grad_phi = gravity.grad_phi(x);
        static_assert(decltype(grad_phi)::size() == 2);

        auto s = cvars_t{0.0,
                         -drho * grad_phi[0],
                         -drho * grad_phi[1],
                         0.0,
                         -(mv[0] * grad_phi[0] + mv[1] * grad_phi[1])};

        return s;
      };

      s += quadrature(s_delta, tri, volume_deg);

      tendency.cvars(i) += s / volume(tri);
    }
  }

  virtual std::string str() const override {
    return "Gravity source term:\n"
           + indent_block(1,
                          string_format("edge degree: %d\nvolume degree: %d\n",
                                        edge_deg,
                                        volume_deg))
      + indent_block(1, type_name<Equilibrium>()) + "\n"
      + indent_block(1, type_name<RC>())
      ;
  }

private:
  euler_t euler;
  std::shared_ptr<Grid> grid;
  std::shared_ptr<grc_t> global_reconstruction;

  int_t edge_deg;
  int_t volume_deg;
};

} // namespace zisa
#endif /* end of include guard */
