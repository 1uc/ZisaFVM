#ifndef REFERENCE_SOLUTION_H_51D7E
#define REFERENCE_SOLUTION_H_51D7E

#include <zisa/grid/grid.hpp>
#include <zisa/grid/point_locator.hpp>
#include <zisa/math/quadrature.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/euler.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {
class ReferenceSolution {
  virtual std::shared_ptr<AllVariables>
  average(const zisa::Grid &coarse_grid) const = 0;
};

template <class Equilibrium>
class EulerReferenceSolution : public ReferenceSolution {
protected:
  using eq_t = Equilibrium;

public:
  EulerReferenceSolution(std::shared_ptr<Grid> fine_grid,
                         std::shared_ptr<AllVariables> fine_vars,
                         const eq_t &eq)
      : fine_grid(std::move(fine_grid)),
        fine_vars(std::move(fine_vars)),
        point_locator(zisa::make_point_locator(fine_grid)),
        n_cvars(fine_vars->cvars.shape(1)) {

    grc = std::make_shared<GlobalReconstruction<Equilibrium, CWENO_AO>>(
        fine_grid, weno_params(), eq);

    grc->compute(*fine_vars);
  }

  virtual std::shared_ptr<AllVariables>
  average(const Grid &coarse_grid) const override {
    LOG_ERR_IF(n_avars == 0, "This needs to be implemented first.")

    auto n_cells = coarse_grid.n_cells;

    auto ref = std::make_shared<AllVariables>(
        AllVariablesDimensions{n_cells, n_cvars, n_avars});

    auto &u_coarse = ref->cvars;

    for (const auto &[i, tri] : triangles(coarse_grid)) {
      for (int_t k_var = 0; k_var < n_cvars; ++k_var) {
        u_coarse(i, k_var) = cvars_average(tri, k_var);
      }
    }

    return ref;
  }

protected:
  double cvars_average(const Triangle &tri, int_t k_var) const {
    auto f = [this, &k_var](const XYZ &x) { return u_ref(x, k_var); };
    return zisa::average(f, tri, quad_deg);
  }

  double u_ref(const XYZ &x, int_t k_var) const {
    auto i_cell = point_locator->locate(x);
    return (*grc)(i_cell)(x)[k_var];
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
  std::shared_ptr<AllVariables> fine_vars;
  std::shared_ptr<PointLocator> point_locator;
  std::shared_ptr<GlobalReconstruction<Equilibrium, CWENO_AO>> grc;

  int_t n_cvars;
  int_t n_avars = 0;

  int_t quad_deg = 4;
};

} // zisa

#endif /* end of include guard */
