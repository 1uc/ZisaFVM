#ifndef EQUILIBRIUM_FLUX_BC_H_YVWNM
#define EQUILIBRIUM_FLUX_BC_H_YVWNM

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class Equilibrium, class EULER>
class EquilibriumFluxBC : public RateOfChange {
private:
  using euler_t = EULER;
  using cvars_t = typename euler_t::cvars_t;

public:
  EquilibriumFluxBC(const euler_t &euler,
                    const Equilibrium &equilibrium,
                    const std::shared_ptr<Grid> &grid,
                    const EdgeRule &qr)
      : euler(euler), equilibrium(equilibrium), grid(grid), qr(qr) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    const auto &eos = euler.eos;

    for (auto &&[e, edge] : exterior_edges(*grid)) {
      auto i = grid->left_right(e).first;
      auto tri = grid->triangle(i);

      auto eq = LocalEquilibrium<Equilibrium>(equilibrium, tri);
      eq.solve(eos.rhoE(cvars_t(current_state.cvars(i))));

      auto flux = [this, &eq, &edge = edge](XY x) {
        auto rhoE = eq.extrapolate(x);

        euler_var_t f = euler.flux(euler.eos.cvars(rhoE));
        inv_coord_transform(f, edge);

        return f;
      };

      auto f = quadrature(qr, flux, edge);
      tendency.cvars(i) -= f / volume(tri);
    }
  }

  virtual std::string str() const override {
    return "Equilibrium flux boundary conditions: \n"
           + indent_block(1, type_name<Equilibrium>().c_str());
  }

private:
  euler_t euler;
  Equilibrium equilibrium;
  std::shared_ptr<Grid> grid;
  EdgeRule qr;
};

} // namespace zisa
#endif /* end of include guard */
