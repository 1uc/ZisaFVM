#ifndef STENCIL_H_FGW9W
#define STENCIL_H_FGW9W

#include <vector>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/math/cone.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/reconstruction/stencil_bias.hpp>
#include <zisa/reconstruction/stencil_params.hpp>

namespace zisa {

class Stencil {
public:
  Stencil() = default;

  /// Simplified constructor for one point stencils.
  explicit Stencil(int_t i_cell);

  Stencil(std::vector<int_t> &l2g,
          const std::shared_ptr<Grid> &grid,
          int_t i_cell,
          const StencilParams &params);

  Stencil(std::vector<int_t> &l2g,
          const std::shared_ptr<Grid> &grid,
          int_t i_cell,
          int_t k,
          const StencilParams &params);

  int_t local(int_t k) const;
  int_t global(int_t k) const;

  /// Factor by which the LSQ problem is over determined.
  double overfit_factor() const;

  /// What is the bias of this stencil (central, one-sided)?
  StencilBias bias() const;

  /// Attainable order with this stencil.
  int order() const;

  /// Requested (maximum) order for this stencil.
  int max_order() const;

  /// Size of the stencil.
  int_t size() const;

protected:
  int_t max_size() const;
  void assign_local_indices(const std::vector<int_t> &global_indices,
                            std::vector<int_t> &l2g);

private:
  int max_order_;
  int order_;

  int_t size_;
  int_t max_size_;

  StencilBias bias_;
  double overfit_factor_;

  array<int_t, 1> local_;
  array<int_t, 1> global_;
};

bool operator==(const Stencil &a, const Stencil &b);
bool operator!=(const Stencil &a, const Stencil &b);

std::vector<int_t> central_stencil(const Grid &grid, int_t i, int_t n_points);

std::vector<int_t>
biased_stencil(const Grid &grid, int_t i_center, int_t k, int_t n_points);

std::vector<int_t>
biased_stencil(const Grid &grid, int_t i, int_t n_points, const Cone &cone);

int deduce_max_order(int_t stencil_size, double factor, int n_dims);
int_t required_stencil_size(int deg, double factor, int n_dims);

} // namespace zisa

#endif /* end of include guard */
