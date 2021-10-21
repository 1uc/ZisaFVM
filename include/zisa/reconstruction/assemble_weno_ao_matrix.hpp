#include <Eigen/Dense>
#include <zisa/grid/grid.hpp>
#include <zisa/reconstruction/stencil.hpp>

namespace zisa {
Eigen::MatrixXd allocate_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

Eigen::MatrixXd
allocate_weno_ao_matrix(const Grid &grid, int order, int n_points);

Eigen::MatrixXd assemble_weno_ao_matrix(const Grid &grid,
                                        const Stencil &stencil);

Eigen::MatrixXd assemble_weno_ao_matrix(
    const Grid &grid, const array_const_view<int_t, 1> &stencil, int order);

void assemble_weno_ao_matrix(Eigen::Ref<Eigen::MatrixXd> &A,
                             const Grid &grid,
                             const Stencil &stencil);

void assemble_weno_ao_matrix(Eigen::Ref<Eigen::MatrixXd> &A,
                             const Grid &grid,
                             const array_const_view<int_t, 1> &stencil,
                             int order);

void assemble_2d_weno_ao_matrix(Eigen::Ref<Eigen::MatrixXd> &A,
                                const Grid &grid,
                                const array_const_view<int_t, 1> &stencil,
                                int order);

void assemble_3d_weno_ao_matrix(Eigen::Ref<Eigen::MatrixXd> &A,
                                const Grid &grid,
                                const array_const_view<int_t, 1> &stencil,
                                int order);
}
