#ifndef ZISA_SELF_GRAVITY_HPP
#define ZISA_SELF_GRAVITY_HPP

#include <zisa/config.hpp>

#include <zisa/grid/grid.hpp>
#include <zisa/math/spherical_shell.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/model/gravity.hpp>
#include <zisa/model/instantaneous_physics.hpp>
#include <zisa/model/poisson_solver.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>

namespace zisa {

class RadialPoissonSolver : public PoissonSolver {
private:
  using euler_t = Euler;

public:
  /// Construct a gravity and poisson solver for it.
  /** NOTE: the `Euler` object does not need to be fully initialized. In
   *   particular, only the equation of state needs to be valid. As part
   *   of constructing the poisson solver, the match gravity will be computed
   *   and set `euler`.
   */
  RadialPoissonSolver(std::shared_ptr<Grid> grid,
                      array<array<int_t, 1>, 1> cell_indices,
                      double gravitational_constant);

  virtual void update(RadialGravity &gravity,
                      const AllVariables &all_variables) const override;

  friend void save(HDF5Writer &writer, const RadialPoissonSolver &solver);

  [[nodiscard]] static RadialPoissonSolver
  load(HDF5Reader &reader, const std::shared_ptr<Grid> &grid);

protected:
  template <class Rho>
  void update_gravity(RadialGravity &gravity, const Rho &rho) const;

  template <class Rho>
  double
  layer_mass(const Rho &rho, const array<double, 1> &radii, int_t layer) const;

private:
  std::shared_ptr<Grid> grid;
  array<array<int_t, 1>, 1> cell_indices;

  double G;
};

void save(HDF5Writer &writer, const RadialPoissonSolver &solver);

///// Solves a radial poisson equation for the gravity.
///** This version does not use the cell-averages directly, instead it uses
/// the
// * (well-balanced) reconstruction to determine the density.
// */
// class CrudeRadialPoissonSolver : public PoissonSolver {
// private:
//  using euler_t = Euler<JankaEOS, RadialGravity>;
//  using cvars_t = euler_t::cvars_t;
//
// public:
//  CrudeRadialPoissonSolver(std::shared_ptr<GlobalReconstruction<cvars_t>>
//  grc,
//                           array<array<int_t, 1>, 1> cell_indices,
//                           array<XYZ, 1> quad_points,
//                           double gravitational_constant)
//      : grc(std::move(grc)),
//        cell_indices(std::move(cell_indices)),
//        quad_points(std::move(quad_points)),
//        G(gravitational_constant) {}
//
//  virtual void update(RadialGravity &gravity,
//                      const AllVariables &all_variables) const override {
//    const auto &radii = gravity.radius_array();
//
//    compute_average_density(radii, all_variables);
//  }
//
// protected:
//  double layer_mass(int_t layer) const;
//  void compute_average_densities(const array<double, 1> &radii) const {
//    int_t n_layers = radii.size() - 1;
//
//    for (int_t l = 0; l < n_layers; ++l) {
//      rho_bar(l) = average_density(radii, l);
//    }
//  }
//
//  double average_density(const array<double, 1> &radii, l) {}
//
// private:
//  std::shared_ptr<GlobalReconstruction<cvars_t>> grc;
//  array<array<int_t, 1>, 1> cell_indices;
//  array<XYZ, 1> quad_points;
//
//  mutable array<double, 1> rho_bar;
//
//  double G;
//};

/// Generate a list of radii such that each layer contains a few cell.
/**
 * @param rel_layer_width width of the layer relative to the local cell-size.
 */
array<double, 1> make_radial_bins(const Grid &grid,
                                  double r_inner,
                                  double r_outer,
                                  double rel_layer_width);

array<double, 1>
make_radial_bins(const Grid &grid, double r_outer, double rel_layer_width);

std::pair<RadialGravity, std::shared_ptr<RadialPoissonSolver>>
make_radial_poisson_solver(const std::shared_ptr<Grid> &grid,
                           double gravitational_constant);

std::shared_ptr<RadialPoissonSolver>
make_radial_poisson_solver(const array<double, 1> &radii,
                           const std::shared_ptr<Grid> &grid,
                           double gravitational_constant);

// std::shared_ptr<CrudeRadialPoissonSolver> make_crude_radial_poisson_solver(
//    const std::shared_ptr<Euler<JankaEOS, RadialGravity>> &euler,
//    const std::shared_ptr<GlobalReconstruction<euler_var_t>> &grc,
//    const Grid &grid,
//    double gravitational_constant,
//    int_t n_quad_points);
}

#endif // ZISA_SELF_GRAVITY_HPP
