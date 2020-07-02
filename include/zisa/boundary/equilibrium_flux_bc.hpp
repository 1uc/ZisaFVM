#ifndef EQUILIBRIUM_FLUX_BC_H_YVWNM
#define EQUILIBRIUM_FLUX_BC_H_YVWNM

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/utils/type_name.hpp>

namespace zisa {

template <class Equilibrium, class EULER>
class EquilibriumFluxBC;

template <class Equilibrium, class EOS, class Gravity>
class EquilibriumFluxBC<Equilibrium, Euler<EOS, Gravity>>
    : public RateOfChange {
private:
  using euler_t = Euler<EOS, Gravity>;
  using cvars_t = typename euler_t::cvars_t;

public:
  EquilibriumFluxBC(std::shared_ptr<euler_t> euler,
                    const Equilibrium &equilibrium,
                    std::shared_ptr<Grid> grid,
                    EdgeRule qr)
      : euler(std::move(euler)),
        equilibrium(equilibrium),
        grid(std::move(grid)),
        qr(std::move(qr)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    const auto &euler = *this->euler;
    const auto &eos = euler.eos;

    for (auto &&[e, face] : exterior_faces(*grid)) {
      auto i = grid->left_right(e).first;
      const auto &cell = grid->cells(i);

      auto eq = LocalEquilibrium<Equilibrium>(equilibrium);
      eq.solve(eos.rhoE(cvars_t(current_state.cvars(i))), cell);

      auto flux = [&euler, &eq, &face = face](XYZ x) {
        auto u = euler.eos.cvars(eq.extrapolate(x));
        auto xvars = euler.eos.xvars(u);
        auto f = euler.flux(u, xvars.p);
        inv_coord_transform(f, face);

        return f;
      };

      auto f = quadrature(face, flux);
      tendency.cvars(i) -= f / volume(cell);
    }
  }

  virtual std::string str() const override {
    return "Equilibrium flux boundary conditions: \n"
           + indent_block(1, type_name<Equilibrium>().c_str());
  }

private:
  std::shared_ptr<euler_t> euler;
  Equilibrium equilibrium;
  std::shared_ptr<Grid> grid;
  EdgeRule qr;
};

} // namespace zisa
#endif /* end of include guard */
