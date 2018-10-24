#include <algorithm>

#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/stencil.hpp>

namespace zisa {

Stencil::Stencil(std::vector<int_t> &l2g,
                 const std::shared_ptr<Grid> &grid,
                 int_t i_cell,
                 const StencilParams &params)
    : max_order_(params.order),
      bias_(deduce_bias(params.bias)),
      overfit_factor_(params.overfit_factor) {

  assert(bias() == StencilBias::central);
  assign_local_indices(central_stencil(*grid, i_cell, max_size()), l2g);

  assert(local_.size() == global_.size());
  order_ = deduce_max_order(local_.size(), overfit_factor());
}

Stencil::Stencil(std::vector<int_t> &l2g,
                 const std::shared_ptr<Grid> &grid,
                 int_t i_cell,
                 int_t k,
                 const StencilParams &params)
    : max_order_(params.order),
      bias_(deduce_bias(params.bias)),
      overfit_factor_(params.overfit_factor) {

  assert(bias() == StencilBias::one_sided);
  assign_local_indices(biased_stencil(*grid, i_cell, k, max_size()), l2g);

  assert(local_.size() == global_.size());
  order_ = deduce_max_order(local_.size(), overfit_factor());
}

void Stencil::assign_local_indices(const std::vector<int_t> &global_indices,
                                   std::vector<int_t> &l2g) {

  assert(global_indices.size() > 0);

  local_ = array<int_t, 1>(shape_t<1>{global_indices.size()});
  global_ = array<int_t, 1>(shape_t<1>{global_indices.size()});

  assert(local_.size() > 0);
  assert(global_.size() > 0);

  for (int_t i = 0; i < global_indices.size(); ++i) {
    auto current = global_indices[i];
    global_(i) = current;

    auto local_index_ptr = std::find(l2g.begin(), l2g.end(), current);
    local_(i) = local_index_ptr - l2g.begin();

    if (local_index_ptr == l2g.end()) {
      l2g.push_back(current);
    }
  }
}

int_t Stencil::local(int_t k) const {
  assert(k < local_.size());
  return local_[k];
}

int_t Stencil::global(int_t k) const {
  assert(k < global_.size());
  return global_[k];
}

int Stencil::order() const { return order_; }
int Stencil::max_order() const { return max_order_; }
StencilBias Stencil::bias() const { return bias_; }
double Stencil::overfit_factor() const { return overfit_factor_; }

int_t Stencil::size() const {
  return required_stencil_size(order() - 1, overfit_factor());
}

int_t Stencil::max_size() const {
  return required_stencil_size(max_order() - 1, overfit_factor());
}

int deduce_max_order(int_t stencil_size, double factor) {
  int deg = 0;
  while (required_stencil_size(deg + 1, factor) <= stencil_size) {
    deg += 1;
  }

  return deg + 1; // order is the degree plus one
}

int_t required_stencil_size(int deg, double factor) {
  if (deg == 0) {
    return 1;
  }

  return int_t(double(poly_dof(deg) - 1) * factor);
}

std::vector<int_t>
central_stencil(const Grid &grid, int_t i_center, int_t n_points) {
  return biased_stencil(
      grid, i_center, n_points, Cone({0.0, 0.0}, {1.0, 0.0}, -1.1));
}

std::vector<int_t>
biased_stencil(const Grid &grid, int_t i_center, int_t k, int_t n_points) {
  const auto &x_center = grid.cell_centers(i_center);
  const auto &v0 = grid.vertex(i_center, k);
  const auto &v1 = grid.vertex(i_center, (k + 1) % grid.max_neighbours);

  return biased_stencil(grid, i_center, n_points, Cone(v0, x_center, v1));
}

std::vector<int_t> biased_stencil(const Grid &grid,
                                  int_t i_center,
                                  int_t n_points,
                                  const Cone &cone) {
  int_t max_points = 5 * n_points;
  int_t max_neighbours = grid.max_neighbours;

  std::vector<int_t> candidates;
  candidates.reserve(3 * max_points);
  candidates.push_back(i_center);

  auto not_found = [&candidates](int_t cand) {
    return std::find(candidates.begin(), candidates.end(), cand)
           == candidates.end();
  };

  auto is_inside = [&cone, &grid](int_t cand) {
    return cone.is_inside(grid.cell_centers(cand));
  };

  for (int_t p = 0; p < max_points; ++p) {
    if (p >= candidates.size()) {
      break;
    }

    int_t j = candidates[p];
    for (int_t k = 0; k < max_neighbours; ++k) {
      if (!grid.is_valid(j, k)) {
        continue;
      }

      int_t cand = grid.neighbours(j, k);

      if (not_found(cand) && is_inside(cand)) {
        candidates.push_back(cand);
      }
    }
  }

  auto closer = [&grid, i_center](int_t i1, int_t i2) {
    const auto &x_center = grid.cell_centers(i_center);
    return zisa::norm(grid.cell_centers(i1) - x_center)
           < zisa::norm(grid.cell_centers(i2) - x_center);
  };

  std::sort(candidates.begin(), candidates.end(), closer);
  candidates.resize(zisa::min(n_points, candidates.size()));
  return candidates;
}

} // namespace zisa
