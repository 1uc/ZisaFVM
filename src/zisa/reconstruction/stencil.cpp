#include <algorithm>

#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/loops/reduction/max.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/reconstruction/stencil.hpp>

namespace zisa {

Stencil::Stencil(int_t i_cell)
    : max_order_(1),
      order_(1),
      size_(1),
      max_size_(1),
      bias_(StencilBias::central),
      overfit_factor_(1.0),
      local_(1),
      global_(1) {
  local_[0] = 0;
  global_[0] = i_cell;
}

Stencil::Stencil(std::vector<int_t> &l2g,
                 const Grid &grid,
                 int_t i_cell,
                 const StencilParams &params)
    : max_order_(params.order),
      bias_(deduce_bias(params.bias)),
      overfit_factor_(params.overfit_factor) {

  int n_dims = grid.n_dims();
  max_size_ = required_stencil_size(max_order() - 1, overfit_factor(), n_dims);

  assert(bias() == StencilBias::central);
  assign_local_indices(central_stencil(grid, i_cell, max_size()), l2g);

  assert(local_.size() == global_.size());
  order_ = deduce_max_order(local_.size(), overfit_factor(), n_dims);
  size_ = required_stencil_size(order() - 1, overfit_factor(), n_dims);
}

Stencil::Stencil(std::vector<int_t> &l2g,
                 const Grid &grid,
                 int_t i_cell,
                 int_t k,
                 const StencilParams &params)
    : max_order_(params.order),
      bias_(deduce_bias(params.bias)),
      overfit_factor_(params.overfit_factor) {
  int n_dims = grid.n_dims();
  max_size_ = required_stencil_size(max_order() - 1, overfit_factor(), n_dims);

  assert(bias() == StencilBias::one_sided);
  assign_local_indices(biased_stencil(grid, i_cell, k, max_size()), l2g);

  assert(local_.size() == global_.size());
  order_ = deduce_max_order(local_.size(), overfit_factor(), n_dims);
  size_ = required_stencil_size(order() - 1, overfit_factor(), n_dims);
}

void Stencil::assign_local_indices(const std::vector<int_t> &global_indices,
                                   std::vector<int_t> &l2g) {

  assert(!global_indices.empty());

  local_ = array<int_t, 1>(shape_t<1>{global_indices.size()});
  global_ = array<int_t, 1>(shape_t<1>{global_indices.size()});

  assert(local_.size() > 0);
  assert(global_.size() > 0);

  for (int_t i = 0; i < global_indices.size(); ++i) {
    auto current = global_indices[i];
    global_(i) = current;

    auto local_index_ptr = std::find(l2g.begin(), l2g.end(), current);
    local_(i) = static_cast<int_t>(local_index_ptr - l2g.begin());

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
int_t Stencil::size() const { return size_; }
int_t Stencil::max_size() const { return max_size_; }

StencilBias Stencil::bias() const { return bias_; }
double Stencil::overfit_factor() const { return overfit_factor_; }
bool operator==(const Stencil &a, const Stencil &b) {
  if (a.size() != b.size()) {
    return false;
  }

  if (a.order() != b.order()) {
    return false;
  }

  if (a.max_order() != b.max_order()) {
    return false;
  }

  if (a.overfit_factor() != b.overfit_factor()) {
    return false;
  }

  if (a.bias() != b.bias()) {
    return false;
  }

  for (int_t i = 0; i < a.size(); ++i) {
    if (a.local(i) != b.local(i)) {
      return false;
    }

    if (a.global(i) != b.global(i)) {
      return false;
    }
  }

  return true;
}

bool operator!=(const Stencil &a, const Stencil &b) { return !(a == b); }

int deduce_max_order(int_t stencil_size, double factor, int n_dims) {
  int deg = 0;
  while (required_stencil_size(deg + 1, factor, n_dims) <= stencil_size) {
    deg += 1;
  }

  return deg + 1; // order is the degree plus one
}

int_t required_stencil_size(int deg, double factor, int n_dims) {
  if (deg == 0) {
    return 1;
  }

  // The polynomials being fitted all have zero mean.
  return int_t(double(poly_dof(deg, n_dims) - 1) * factor);
}

std::vector<int_t>
central_stencil(const Grid &grid, int_t i_center, int_t n_points) {
  return biased_stencil(
      grid, i_center, n_points, Cone({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, -1.1));
}

static double compute_cos_alpha(const Grid &grid, int_t i, int_t k) {
  auto max_neighbours = grid.max_neighbours;
  auto element_type = (max_neighbours == 3 ? GMSHElementType::triangle
                                           : GMSHElementType::tetrahedron);

  const auto &x_center = grid.cell_centers(i);

  auto range = index_range(GMSHElementInfo::n_vertices_per_face(element_type));

  double r_max = zisa::reduce::max(
      range, [&grid, &x_center, element_type, i, k](int_t rel) {
        auto k_rel
            = GMSHElementInfo::relative_vertex_index(element_type, k, rel);
        const auto &vi = grid.vertices(grid.vertex_indices(i, k_rel));

        return zisa::norm(x_center - vi);
      });

  double r = zisa::norm(grid.face_center(i, k) - x_center);

  return r / r_max;
}

std::vector<int_t>
biased_stencil(const Grid &grid, int_t i_center, int_t k, int_t n_points) {
  const auto &cell_center = grid.cell_centers(i_center);
  const auto &face_center = grid.face_center(i_center, k);

  double cos_alpha = compute_cos_alpha(grid, i_center, k);

  auto v = XYZ(face_center - cell_center);
  auto v_hat = XYZ(normalize(v));

  return biased_stencil(
      grid, i_center, n_points, Cone(cell_center, v_hat, cos_alpha));
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
