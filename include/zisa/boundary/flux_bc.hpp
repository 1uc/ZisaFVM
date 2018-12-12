#ifndef FLUX_BC_H_KDUG6
#define FLUX_BC_H_KDUG6

#include <zisa/grid/grid.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

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
                       const AllVariables &current_state,
                       double /* t */) const override {

    for (auto &&[e, edge] : exterior_edges(*grid)) {
      auto i = grid->left_right(e).first;

      auto u = cvars_t(current_state.cvars(i));
      coord_transform(u, edge);

      auto f = model.flux(u);
      inv_coord_transform(f, edge);

      tendency.cvars(i) -= volume(edge) / grid->volumes(i) * f;
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
