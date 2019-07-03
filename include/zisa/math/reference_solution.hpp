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
};

template <class Equilibrium, class Scaling>
class EulerReferenceSolution : public ReferenceSolution {
protected:
  using eq_t = Equilibrium;
  using grc_t = EulerGlobalReconstruction<Equilibrium, CWENO_AO, Scaling>;

public:
  EulerReferenceSolution(std::shared_ptr<Grid> fine_grid,
                         const std::shared_ptr<AllVariables> &fine_vars,
                         const eq_t &eq,
                         const Scaling &scaling)
      : fine_grid(std::move(fine_grid)), n_cvars(fine_vars->cvars.shape(1)) {

    LOG_ERR_IF(fine_vars->dims().n_avars != n_avars, "Dimension mismatch.");

    grc = std::make_shared<grc_t>(this->fine_grid, weno_params(), eq, scaling);
    grc->compute(*fine_vars);
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
        const auto &cell = coarse_grid.cells(i);
        auto u_bar = cvars_average(cell, i_guess);
        for (int_t k_var = 0; k_var < n_cvars; ++k_var) {
          u_coarse(i, k_var) = u_bar[k_var];
        }
      }
    }

    return ref;
  }

protected:
  euler_var_t cvars_average(const Cell &cell, int_t &i_guess) const {
    auto f = [this, &i_guess](const XYZ &x) { return u_ref(x, i_guess); };
    return zisa::average(cell, f);
  }

  euler_var_t u_ref(const XYZ &x, int_t &i_guess) const {
    // auto i_cell = point_locator->locate(x);
    int_t max_iter = fine_grid->n_cells;
    auto i_cell = locate(*fine_grid, x, i_guess, max_iter);

    LOG_ERR_IF(!i_cell, "Failed to locate the cell.");
    i_guess = *i_cell;

    return euler_var_t((*grc)(*i_cell)(x));
  }

  HybridWENOParams weno_params() {
    return HybridWENOParams{
        {{5, 2, 2, 2}, {"c", "b", "b", "b"}, {2.0, 1.5, 1.5, 1.5}},
        {100.0, 1.0, 1.0, 1.0},
        1e-10,
        4.0};
  }

private:
  std::shared_ptr<Grid> fine_grid;
  std::shared_ptr<grc_t> grc;

  int_t n_cvars;
  int_t n_avars = 0;

  int_t quad_deg = 4;

  mutable int_t last_index_ = 0;
};

} // zisa

#endif /* end of include guard */
