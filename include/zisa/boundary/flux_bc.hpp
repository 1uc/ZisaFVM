#ifndef FLUX_BC_H_KDUG6
#define FLUX_BC_H_KDUG6

#include <zisa/grid/grid.hpp>
#include <zisa/ode/rate_of_change.hpp>

namespace zisa {

template <class Model, class RC>
class FluxBC : public RateOfChange {
private:
  using cvars_t = typename Model::cvars_t;

public:
  FluxBC(const Model &model,
         const EdgeRule &edge_rule,
         const std::shared_ptr<Grid> &grid,
         const std::shared_ptr<GlobalReconstruction<RC>> &global_reconstruction)
      : model(model),
        edge_rule(edge_rule),
        grid(grid),
        global_reconstruction(global_reconstruction) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables & /* current_state */,
                       double /* t */) const override {

    for (auto &&[e, edge] : exterior_edges(*grid)) {
      auto i = grid->left_right(e).first;

      auto flux = [this, i, &edge](const XY &x) -> cvars_t {
        cvars_t u;
        for (int_t k = 0; k < cvars_t::size(); ++k) {
          u[k] = (*global_reconstruction)(i, k)(x);
        }
        coord_transform(u, edge);

        return model.flux(u);
      };

      auto f = quadrature(edge_rule, flux, edge);
      inv_coord_transform(f, edge);

      tendency.cvars(i) -= f / grid->volumes(i);
    }
  }

  virtual std::string str() const override {
    return "Flux boundary conditions.";
  }

private:
  Model model;
  EdgeRule edge_rule;

  std::shared_ptr<Grid> grid;
  std::shared_ptr<GlobalReconstruction<RC>> global_reconstruction;
};

} // namespace zisa

#endif /* end of include guard */
