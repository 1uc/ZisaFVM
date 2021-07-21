// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef FLUX_LOOP_H_BWHPN
#define FLUX_LOOP_H_BWHPN

#include <zisa/grid/grid.hpp>
#include <zisa/math/edge_rule.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/local_eos_state.hpp>
#include <zisa/ode/rate_of_change.hpp>
#include <zisa/parallelization/halo_exchange.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>
#include <zisa/utils/indent_block.hpp>
#include <zisa/utils/timer.hpp>

namespace zisa {

template <class Model, class Flux, class LEOS, class GRC>
class FluxLoop : public RateOfChange {
protected:
  using cvars_t = typename Model::cvars_t;
  using eos_t = typename LEOS::eos_t;
  using grc_t = GRC;

public:
  FluxLoop(std::shared_ptr<Grid> grid_,
           std::shared_ptr<Model> model_,
           std::shared_ptr<LEOS> local_eos_,
           std::shared_ptr<grc_t> global_reconstruction_,
           std::shared_ptr<HaloExchange> halo_exchange_,
           EdgeRule edge_rule)
      : grid(std::move(grid_)),
        model(std::move(model_)),
        local_eos(std::move(local_eos_)),
        global_reconstruction(std::move(global_reconstruction_)),
        edge_rule(std::move(edge_rule)),
        halo_exchange(std::move(halo_exchange_)) {

    avar_flux_allocator
        = std::make_shared<block_allocator<array<double, 1>>>(128);

    interior_cells.reserve(this->grid->n_cells);
    exterior_cells.reserve(this->grid->n_cells / 2);

    for (int_t i = 0; i < this->grid->n_cells; ++i) {
      const auto &s = this->global_reconstruction->stencil(i);
      for (int_t j : s) {
        if (this->grid->cell_flags[j].ghost_cell) {
          exterior_cells.push_back(i);
          break;
        }
      }

      if (exterior_cells.empty() || exterior_cells.back() != i) {
        interior_cells.push_back(i);
      }
    }

    std::sort(interior_cells.begin(), interior_cells.end());
    std::sort(exterior_cells.begin(), exterior_cells.end());

    interior_faces.reserve(this->grid->n_interior_edges);
    exterior_faces.reserve((this->grid->n_edges - this->grid->n_interior_edges)
                           * 2);

    // In the context of a grid, any edge which has two neighbours is an
    // interior edge. Here we attempt to find the edges which are safe to update
    // before the halo has been exchanged.
    for (int_t e = 0; e < this->grid->n_interior_edges; ++e) {
      auto [iL, iR] = this->grid->left_right(e);

      auto is_interior
          = !contains(exterior_cells, iL) && !contains(exterior_cells, iR);

      if (is_interior) {
        interior_faces.push_back(e);

      } else {
        auto is_left_ghost = this->grid->cell_flags[iL].ghost_cell;
        auto is_right_ghost = this->grid->cell_flags[iR].ghost_cell;

        if (!(is_left_ghost && is_right_ghost)) {
          exterior_faces.push_back(e);
        }
      }
    }
  }

  bool contains(const std::vector<int_t> &cells, int_t i) const {
    return std::lower_bound(cells.begin(), cells.end(), i) != cells.end();
  }

  virtual void compute(AllVariables &tendency,
                       const AllVariables &current_state,
                       double /* t */) const override {

    (*halo_exchange)(const_cast<AllVariables &>(current_state));
    compute_patch(tendency, current_state, interior_cells, interior_faces);
    (*halo_exchange).wait();
    compute_patch(tendency, current_state, exterior_cells, exterior_faces);
  }

  void compute_patch(AllVariables &tendency,
                     const AllVariables &current_state,
                     const array_const_view<int_t, 1> &cells,
                     const array_const_view<int_t, 1> &edges) const {

    (*local_eos).compute(current_state, cells);
    (*global_reconstruction).compute(current_state, cells);

    const auto n_avars = tendency.avars.shape(1);
    const auto n_interior_edges = edges.size();
#if ZISA_HAS_OPENMP == 1
#pragma omp parallel
#endif
    {
      auto qnf_guard = fetch_avars_flux_buffer(n_avars);
      auto &qnf = *qnf_guard;

#if ZISA_HAS_OPENMP == 1
#pragma omp for ZISA_OMP_FOR_SCHEDULE_DEFAULT
#endif
      for (int_t ie = 0; ie < n_interior_edges; ++ie) {
        auto e = edges[ie];

        const auto &face = grid->faces(e);

        int_t iL, iR;
        std::tie(iL, iR) = grid->left_right(e);
        const auto &eosL = *(*local_eos)(iL);
        const auto &eosR = *(*local_eos)(iR);

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

          const auto [f, speeds] = numerical_flux(eosL, uL, eosR, uR);
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
  auto numerical_flux(const eos_t &eosL,
                      const cvars_t &uL,
                      const eos_t &eosR,
                      const cvars_t &uR) const {
    return Flux::flux(*model, eosL, uL, eosR, uR);
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

  std::shared_ptr<LEOS> local_eos;
  std::shared_ptr<grc_t> global_reconstruction;
  mutable std::shared_ptr<block_allocator<array<double, 1>>>
      avar_flux_allocator;

  EdgeRule edge_rule;

  std::shared_ptr<HaloExchange> halo_exchange;

  std::vector<int_t> interior_cells;
  std::vector<int_t> exterior_cells;
  std::vector<int_t> interior_faces;
  std::vector<int_t> exterior_faces;
};

} // namespace zisa

#endif
