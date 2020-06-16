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
        edge_rule(std::move(edge_rule)) {

    avar_flux_allocator
        = std::make_shared<block_allocator<array<double, 1>>>(128);
  }

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    (*global_reconstruction).compute(current_state);

    const auto n_avars = tendency.avars.shape(1);
    const auto n_interior_edges = grid->n_interior_edges;
#if ZISA_HAS_OPENMP == 1
#pragma omp parallel
#endif
    {
      auto qnf_guard = fetch_avars_flux_buffer(n_avars);
      auto &qnf = *qnf_guard;

#if ZISA_HAS_OPENMP == 1
#pragma omp for ZISA_OMP_FOR_SCHEDULE_DEFAULT
#endif
      for (int_t e = 0; e < n_interior_edges; ++e) {
        const auto &face = grid->faces(e);

        int_t iL, iR;
        std::tie(iL, iR) = grid->left_right(e);

        auto rc = [this, &face = face](int_t i, const XYZ &x) {
          auto u = cvars_t((*global_reconstruction)(i)(x));
          coord_transform(u, face);

          return u;
        };

        auto nf = cvars_t::zeros();
        zisa::fill(qnf, 0.0);

        const auto n_qr = face.qr.weights.size();
        for (int_t k = 0; k < n_qr; ++k) {
          const auto w = face.qr.weights[k];
          const auto x = face.qr.points[k];

          const auto uL = rc(iL, x);
          const auto uR = rc(iR, x);

          const auto [f, speeds] = numerical_flux(uL, uR);
          nf += w * f;

          for (int_t a = 0; a < n_avars; ++a) {
            const auto qL = (*global_reconstruction)(iL).tracer(x, a);
            const auto qR = (*global_reconstruction)(iR).tracer(x, a);
            qnf[a] += w * tracer_flux(uL, uR, qL, qR, speeds);
          }
        }

        inv_coord_transform(nf, face);

        for (int_t k = 0; k < cvars_t::size(); ++k) {
          auto nfL = nf(k) / grid->volumes(iL);
#if ZISA_HAS_OPENMP == 1
#pragma omp atomic
#endif
          tendency.cvars(iL, k) -= nfL;

          auto nfR = nf(k) / grid->volumes(iR);
#if ZISA_HAS_OPENMP == 1
#pragma omp atomic
#endif
          tendency.cvars(iR, k) += nfR;
        }

        for (int_t k = 0; k < n_avars; ++k) {
          const auto qfL = qnf(k) / grid->volumes(iL);
#if ZISA_HAS_OPENMP == 1
#pragma omp atomic
#endif
          tendency.avars(iL, k) -= qfL;

          const auto qfR = qnf(k) / grid->volumes(iR);
#if ZISA_HAS_OPENMP == 1
#pragma omp atomic
#endif
          tendency.avars(iR, k) += qfR;
        }
      }
    }
  }

  locked_ptr<array<double, 1>> fetch_avars_flux_buffer(int_t n_avars) const {
    return avar_flux_allocator->allocate(shape_t<1>(n_avars));
  }

  virtual std::string str() const override {
    auto block = model->str() + "\n" + edge_rule.str() + "\n"
                 + global_reconstruction->str();

    return "Flux loop: \n" + indent_block(1, block);
  }

private:
  auto numerical_flux(const cvars_t &uL, const cvars_t &uR) const {
    return Flux::flux(*model, uL, uR);
  }

  auto tracer_flux(const cvars_t &uL,
                   const cvars_t &uR,
                   const double qL,
                   const double qR,
                   const std::tuple<double, double, double> &speeds) const {
    return Flux::tracer_flux(uL, uR, qL, qR, speeds);
  }

private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<Model> model;

  std::shared_ptr<grc_t> global_reconstruction;
  mutable std::shared_ptr<block_allocator<array<double, 1>>>
      avar_flux_allocator;

  EdgeRule edge_rule;
};

} // namespace zisa

#endif
