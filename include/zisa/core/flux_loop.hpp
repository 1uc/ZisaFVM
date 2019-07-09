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

template <class Model, class Flux, class GRC>
class FluxLoop : public RateOfChange {
protected:
  using cvars_t = typename Model::cvars_t;
  using grc_t = GRC;

public:
  FluxLoop(std::shared_ptr<Grid> grid,
           std::shared_ptr<Model> model,
           std::shared_ptr<grc_t> global_reconstruction,
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

#pragma omp parallel for ZISA_OMP_FOR_SCHEDULE_DEFAULT
    for (int_t e = 0; e < n_interior_edges; ++e) {
      auto face = grid->faces(e);

      int_t iL, iR;
      std::tie(iL, iR) = grid->left_right(e);

      auto rc = [this, &face = face](int_t i, const XYZ &x) {
        auto u = cvars_t((*global_reconstruction)(i)(x));
        coord_transform(u, face);

        return u;
      };

      auto flux = [this, &rc, iL = iL, iR = iR](const XYZ &x) -> cvars_t {
        auto uL = rc(iL, x);
        auto uR = rc(iR, x);

        return numerical_flux(uL, uR);
      };

      auto nf = quadrature(face, flux);
      inv_coord_transform(nf, face);

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
    auto block = model->str() + "\n" + edge_rule.str() + "\n"
                 + global_reconstruction->str();

    return "Flux loop: \n" + indent_block(1, block);
  }

private:
  cvars_t numerical_flux(const cvars_t &uL, const cvars_t &uR) const {
    return Flux::flux(*model, uL, uR);
  }

private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<Model> model;

  std::shared_ptr<grc_t> global_reconstruction;

  EdgeRule edge_rule;
};

} // namespace zisa

#endif
