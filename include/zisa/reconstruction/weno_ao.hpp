/*
 *
 */

#ifndef WENO_AO_H_A1UM2
#define WENO_AO_H_A1UM2

#include <Eigen/Dense>

#include <zisa/config.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/cone.hpp>
#include <zisa/memory/array_view.hpp>

namespace zisa {

struct LocalIndex {
  int_t local;
  int_t global;
};

using LocalBuffer = array<double, 1>;

class WENO_AO {
public:
  static constexpr int_t max_degree() { return 4ul; }

public:
  using QR = Eigen::FullPivHouseholderQR<Eigen::MatrixXd>;

public:
  WENO_AO(const Grid &grid, int_t i_cell);
  Poly2D<4> reconstruct(const LocalBuffer &buffer) const;

protected:
  array<LocalIndex, 1>
  assign_local_indices(const std::vector<int_t> &global_indices,
                       std::vector<int_t> &l2g);

  void compute_stencils(const Grid &grid, int_t i_cell);
  void compute_qr(const Grid &grid, int_t i_cells);

private:
  array<array<LocalIndex, 1>, 1> stencils;
  array<QR, 1> qr;
  int order;
};

int_t deduce_max_order(const std::vector<int_t> &stencil, double factor);
int_t required_stencil_size(int_t deg, double factor);

array<double, 1>
normalized_moments(const Triangle &tri, double length, int_t deg);

std::vector<int_t> central_stencil(const Grid &grid, int_t i, int_t n_points);

std::vector<int_t>
biased_stencil(const Grid &grid, int_t i_center, int_t k, int_t n_points);

std::vector<int_t>
biased_stencil(const Grid &grid, int_t i, int_t n_points, const Cone &cone);

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const array<LocalIndex, 1> &stencil,
                                        int_t order,
                                        double factor);

} // namespace zisa
#endif /* end of include guard */
