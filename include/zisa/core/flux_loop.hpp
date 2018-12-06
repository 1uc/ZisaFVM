#ifndef FLUX_LOOP_H_BWHPN
#define FLUX_LOOP_H_BWHPN

#include <zisa/grid/grid.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>
#include <zisa/utils/indent_block.hpp>

namespace zisa {

template <class RC, class Model, class Flux>
class FluxLoop : public RateOfChange {
protected:
  using cvars_t = typename Model::cvars_t;

public:
  FluxLoop(std::shared_ptr<Grid> grid,
           Model model,
           std::shared_ptr<GlobalReconstruction<RC>> global_reconstruction,
           EdgeRule edge_rule)
      : grid(std::move(grid)),
        model(std::move(model)),
        global_reconstruction(std::move(global_reconstruction)),
        edge_rule(std::move(edge_rule)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    (*global_reconstruction).compute(current_state);

    for (const auto &[e, edge] : interior_edges(*grid)) {
      auto [iL, iR] = grid->left_right(e);

      auto rc = [this, &edge = edge](int_t i, const XY &x) {
        cvars_t u;
        for (int_t k = 0; k < cvars_t::size(); ++k) {
          u[k] = (*global_reconstruction)(i, k)(x);
        }
        coord_transform(u, edge);

        return u;
      };

      auto flux = [this, &rc, iL = iL, iR = iR](const XY &x) -> cvars_t {
        auto uL = rc(iL, x);
        auto uR = rc(iR, x);

        return numerical_flux(uL, uR);
      };

      auto nf = quadrature(edge_rule, flux, edge);
      inv_coord_transform(nf, edge);

      tendency.cvars(iL) -= nf / grid->volumes(iL);
      tendency.cvars(iR) += nf / grid->volumes(iR);
    }
  }

  virtual std::string str() const override {
    auto block = model.str() + "\n" + edge_rule.str() + "\n"
                 + global_reconstruction->str();

    return "Flux loop: \n" + indent_block(1, block);
  }

private:
  cvars_t numerical_flux(const cvars_t &uL, const cvars_t &uR) const {
    return Flux::flux(model, uL, uR);
  }

private:
  std::shared_ptr<Grid> grid;
  Model model;
  std::shared_ptr<GlobalReconstruction<RC>> global_reconstruction;

  EdgeRule edge_rule;
};

} // namespace zisa

#endif
