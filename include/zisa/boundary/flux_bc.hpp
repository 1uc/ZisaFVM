#ifndef FLUX_BC_H_KDUG6
#define FLUX_BC_H_KDUG6

#include <zisa/grid/grid.hpp>
#include <zisa/ode/rate_of_change.hpp>

namespace zisa {

template <class Model>
class FluxBC : public RateOfChange {
private:
  using cvars_t = typename Model::cvars_t;

public:
  FluxBC(const Model &model, const std::shared_ptr<Grid> &grid)
      : model(model), grid(grid) {}

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
  std::shared_ptr<Grid> grid;
};

} // namespace zisa

#endif /* end of include guard */
