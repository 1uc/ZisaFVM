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

    auto n_interior_edges = grid->n_interior_edges;

#pragma omp parallel for schedule(guided)
    for (int_t e = 0; e < n_interior_edges; ++e) {
      auto edge = grid->edge(e);

      int_t iL, iR;
      std::tie(iL, iR) = grid->left_right(e);

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

      for (int_t k = 0; k < cvars_t::size(); ++k) {
        auto nfL = nf(k) / grid->volumes(iL);
#pragma omp atomic
        tendency.cvars(iL, k) -= nfL;

        auto nfR = nf(k) / grid->volumes(iR);
#pragma omp atomic
        tendency.cvars(iR, k) += nfR;
      }
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
