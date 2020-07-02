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
  FluxBC(std::shared_ptr<Model> model, std::shared_ptr<Grid> grid)
      : model(std::move(model)), grid(std::move(grid)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    for (auto &&[e, face] : exterior_faces(*grid)) {
      auto i = grid->left_right(e).first;

      auto u = cvars_t(current_state.cvars(i));
      coord_transform(u, face);

      auto xvars = model->eos.xvars(u);
      auto f = model->flux(u, xvars.p);
      inv_coord_transform(f, face);

      tendency.cvars(i) -= volume(face) / grid->volumes(i) * f;
    }
  }

  virtual std::string str() const override {
    return "Flux boundary conditions.";
  }

private:
  std::shared_ptr<Model> model;
  std::shared_ptr<Grid> grid;
};

} // namespace zisa

#endif /* end of include guard */
