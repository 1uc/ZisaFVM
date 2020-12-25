#ifndef FLUX_BC_H_KDUG6
#define FLUX_BC_H_KDUG6

#include <zisa/grid/grid.hpp>
#include <zisa/ode/rate_of_change.hpp>

namespace zisa {

template <class EOS>
class FluxBC : public RateOfChange {
private:
  using cvars_t = typename Euler::cvars_t;

public:
  FluxBC(std::shared_ptr<Euler> euler,
         std::shared_ptr<LocalEOSState<EOS>> local_eos,
         std::shared_ptr<Grid> grid)
      : euler(std::move(euler)),
        local_eos(std::move(local_eos)),
        grid(std::move(grid)) {}

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    for (auto &&[e, face] : exterior_faces(*grid)) {
      auto i = grid->left_right(e).first;
      const auto &eos = (*local_eos)(i);

      auto u = cvars_t(current_state.cvars(i));
      coord_transform(u, face);

      auto xvars = eos->xvars(u);
      auto f = euler->flux(u, xvars.p);
      inv_coord_transform(f, face);

      tendency.cvars(i) -= volume(face) / grid->volumes(i) * f;
    }
  }

  virtual std::string str() const override {
    return "Flux boundary conditions.";
  }

private:
  std::shared_ptr<Euler> euler;
  std::shared_ptr<LocalEOSState<EOS>> local_eos;
  std::shared_ptr<Grid> grid;
};


class NoFluxBC : public RateOfChange {
public:
  virtual void compute(AllVariables &,
                       const AllVariables &,
                       double /* t */) const override {
    return;
  }

  virtual std::string str() const override {
    return "No flux boundary conditions.";
  }

};

} // namespace zisa

#endif /* end of include guard */
