#ifndef REFERENCE_SOLUTION_H_51D7E
#define REFERENCE_SOLUTION_H_51D7E

#include <zisa/grid/grid.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {
class ReferenceSolution {
public:
  virtual ~ReferenceSolution() = default;

  virtual std::shared_ptr<AllVariables>
  average(const zisa::Grid &coarse_grid) const = 0;

  virtual double q_ref(const XYZ &x, int_t k, int_t &i_guess) const = 0;
};

template <class Equilibrium, class Scaling>
class EulerReferenceSolution : public ReferenceSolution {
protected:
  using eq_t = Equilibrium;
  using grc_t = EulerGlobalReconstruction<Equilibrium, CWENO_AO, Scaling>;

public:
  EulerReferenceSolution(std::shared_ptr<Grid> fine_grid,
                         const AllVariables &fine_vars,
                         std::shared_ptr<grc_t> grc)
      : fine_grid(std::move(fine_grid)),
        grc(std::move(grc)),
        n_cvars(fine_vars.cvars.shape(1)) {

    // FIXME we now have tracers.
    LOG_ERR_IF(fine_vars.dims().n_avars != n_avars, "Dimension mismatch.");
    grc->compute(fine_vars);
  }

  virtual std::shared_ptr<AllVariables>
  average(const Grid &coarse_grid) const override {
    auto n_cells = coarse_grid.n_cells;

    auto ref = std::make_shared<AllVariables>(
        AllVariablesDimensions{n_cells, n_cvars, n_avars});

    auto &u_coarse = ref->cvars;

    {
      int_t i_guess = 0;

      for (int_t i = 0; i < n_cells; ++i) {
        if (coarse_grid.cell_flags(i).ghost_cell) {
          for (int_t k_var = 0; k_var < n_cvars; ++k_var) {
            u_coarse(i, k_var) = 0.0;
          }
        } else {
          const auto &cell = coarse_grid.cells(i);
          auto u_bar = cvars_average(cell, i_guess);
          for (int_t k_var = 0; k_var < n_cvars; ++k_var) {
            u_coarse(i, k_var) = u_bar[k_var];
          }
        }
      }
    }

    return ref;
  }

  double q_ref(const XYZ &x, int_t k, int_t &i_guess) const override {
    auto i_cell = locate(*fine_grid, x, i_guess);
    LOG_ERR_IF(!i_cell,
               string_format("Failed to locate the cell. x = %s",
                             format_as_list(x).c_str()));

    i_guess = *i_cell;
    return (*grc)(*i_cell)(x)[k];
  }

protected:
  euler_var_t cvars_average(const Cell &cell, int_t &i_guess) const {
    auto f = [this, &i_guess](const XYZ &x) { return u_ref(x, i_guess); };
    return zisa::average(cell, f);
  }

  euler_var_t u_ref(const XYZ &x, int_t &i_guess) const {
    auto i_cell = locate(*fine_grid, x, i_guess);
    LOG_ERR_IF(!i_cell,
               string_format("Failed to locate the cell. x = %s",
                             format_as_list(x).c_str()));

    i_guess = *i_cell;

    return euler_var_t((*grc)(*i_cell)(x));
  }

  HybridWENOParams weno_params() {
    if (fine_grid->n_dims() == 2) {
      return HybridWENOParams{
          {{5, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
          {100.0, 1.0, 1.0, 1.0},
          1e-10,
          4.0};
    }

    if (fine_grid->n_dims() == 3) {
      return HybridWENOParams{{{4, 2, 2, 2, 2},
                               {"c", "b", "b", "b", "b"},
                               {3.0, 2.0, 2.0, 2.0, 2.0}},
                              {100.0, 1.0, 1.0, 1.0, 1.0},
                              1e-10,
                              4.0};
    }

    LOG_ERR("Implement first.");
  }

private:
  std::shared_ptr<Grid> fine_grid;
  std::shared_ptr<grc_t> grc;

  int_t n_cvars;
  int_t n_avars = 0;

  mutable int_t last_index_ = 0;
};

} // zisa

#endif /* end of include guard */
