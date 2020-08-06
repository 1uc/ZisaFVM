#include <algorithm>

#include <random>
#include <zisa/grid/gmsh_reader.hpp>
#include <zisa/loops/reduction/max.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/math/tetrahedral_rule.hpp>
#include <zisa/math/triangular_rule.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/reconstruction/lsq_solver.hpp>
#include <zisa/reconstruction/stencil.hpp>
#include <zisa/utils/indent_block.hpp>

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

Stencil::Stencil(std::vector<int_t> &global_indices,
                 const StencilParams &params)
    : max_order_(params.order),
      bias_(deduce_bias(params.bias)),
      overfit_factor_(params.overfit_factor),
      local_(0),
      global_(global_indices.size()) {

  global_ = array_const_view(global_indices);
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
  assign_local_indices(biased_stencil(grid, i_cell, k, max_size(), max_order_),
                       l2g);

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
  return int_t(double(poly_dof(deg, n_dims) - 1) * factor + 1);
}

// --- Code for generating stencils. -------------------------------------------
DenormalizedRule make_stencil_selection_query_points(const Grid &grid,
                                                     int_t i) {

  if (grid.is_triangular()) {
    auto qr_hat = cached_triangular_quadrature_rule(MAX_TRIANGULAR_RULE_DEGREE);
    return denormalize(qr_hat, triangle(grid, i));
  } else if (grid.is_tetrahedral()) {
    auto qr_hat = cached_tetrahedral_rule(MAX_TETRAHEDRAL_RULE_DEGREE);
    return denormalize(qr_hat, tetrahedron(grid, i));
  }

  LOG_ERR("Implement the missing case.");
}

std::vector<int_t> region_based_candidates(const Grid &grid,
                                           int_t i_center,
                                           int_t n_points,
                                           const Region &region) {

  int_t max_points = 5 * n_points;
  int_t max_neighbours = grid.max_neighbours;

  std::vector<int_t> candidates;
  candidates.reserve(3 * max_points);
  candidates.push_back(i_center);

  auto not_found = [&candidates](int_t cand) {
    return std::find(candidates.begin(), candidates.end(), cand)
           == candidates.end();
  };

  auto is_inside = [&region, &grid](int_t cand) {
    auto qr = make_stencil_selection_query_points(grid, cand);
    for (const auto &x : qr.points) {
      if (region.is_inside(x)) {
        return true;
      }
    }

    return region.is_inside(grid.cell_centers(cand));
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

  return candidates;
}

std::vector<int_t> region_based_stencil(const Grid &grid,
                                        int_t i_center,
                                        int_t n_points,
                                        const Region &region) {

  auto candidates = region_based_candidates(grid, i_center, n_points, region);

  auto closer = [&grid, i_center](int_t i1, int_t i2) {
    const auto &x_center = grid.cell_centers(i_center);
    return zisa::norm(grid.cell_centers(i1) - x_center)
           < zisa::norm(grid.cell_centers(i2) - x_center);
  };

  std::sort(candidates.begin(), candidates.end(), closer);
  candidates.resize(zisa::min(n_points, candidates.size()));
  return candidates;
}

std::shared_ptr<Region>
make_cone(const Grid &grid, int_t i_center, const XYZ &cone_center, int_t k) {

  auto element_type = grid.element_type();

  auto rel_vertex = [&grid, &element_type](int_t i, int_t k, int_t rel) {
    auto kk = GMSHElementInfo::relative_vertex_index(element_type, k, rel);
    return grid.vertices[grid.vertex_indices(i, kk)];
  };

  if (element_type == GMSHElementType::triangle) {
    return std::make_shared<TriangularCone>(
        cone_center, rel_vertex(i_center, k, 0), rel_vertex(i_center, k, 1));
  } else {
    return std::make_shared<TetrahedralCone>(cone_center,
                                             rel_vertex(i_center, k, 0),
                                             rel_vertex(i_center, k, 1),
                                             rel_vertex(i_center, k, 2));
  }
}

std::vector<int_t> conservative_stencil(const Grid &grid,
                                        int_t i_center,
                                        int_t k,
                                        int_t n_points) {

  auto element_type = grid.element_type();
  const auto &off_vertex = grid.vertices(grid.vertex_indices(
      i_center, GMSHElementInfo::relative_off_vertex_index(element_type, k)));

  auto region = make_cone(grid, i_center, off_vertex, k);
  return region_based_stencil(grid, i_center, n_points, *region);
}

std::vector<int_t> less_conservative_stencil(const Grid &grid,
                                             int_t i_center,
                                             int_t k,
                                             int_t n_points) {

  const auto &cell_center = grid.cell_centers(i_center);
  auto region = make_cone(grid, i_center, cell_center, k);
  return region_based_stencil(grid, i_center, n_points, *region);
}

std::vector<int_t> tryhard_stencil(
    const Grid &grid,
    const std::function<bool(const array_const_view<int_t, 1> &)> &is_good,
    int_t i_center,
    int_t k,
    int_t n_points) {

  // Idea:
  // Everything right of face k is a candidate, keep picking
  // random stencils until something works.

  auto cone_center = grid.face_center(i_center, k);
  auto region = make_cone(grid, i_center, cone_center, k);

  auto candidates = region_based_candidates(grid, i_center, n_points, *region);
  if (candidates.size() < n_points) {
    return std::vector<int_t>{i_center};
  }
  if (!is_good(candidates)) {
    std::stringstream ss;
    ss << "It is impossible to find a suitable stencil.\n";
    ss << "  candidates = " << format_as_list(candidates) << "\n";
    ss << "  i_center = " << i_center << "\n";
    ss << "  k = " << k << "\n";
    ss << "  n_points = " << n_points << "\n";
    LOG_ERR(ss.str());
  }

  std::random_device rd;
  for (int iter = 0; iter < 100; ++iter) {
    std::shuffle(candidates.begin() + 1, candidates.end(), std::mt19937(rd()));

    if (is_good(array_const_view(shape_t<1>(n_points), candidates.data()))) {
      candidates.resize(n_points);
      LOG_WARN_IF(
          iter >= 1,
          string_format("It took %d iterations to find a stencil.", iter + 1));
      return candidates;
    }
  }
  LOG_ERR("Did not find a suitable stencil given the candidates.");
}

std::vector<int_t> biased_stencil(
    const Grid &grid, int_t i_center, int_t k, int_t n_points, int order) {

  auto is_good = [&](const array_const_view<int_t, 1> &s) {
    auto A = assemble_weno_ao_matrix(grid, s, order);
    auto svd = Eigen::JacobiSVD(A);

    if (svd.rank() != A.cols()) {

      Eigen::IOFormat python_fmt(
          Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");

      std::stringstream ss;
      ss.precision(16);
      ss << "n_dims = " << grid.n_dims() << "\n";
      ss << "rank = " << svd.rank() << "\n";
      ss << "cols = " << svd.cols() << "\n";
      ss << "A = \n" << A.format(python_fmt) << "\n";
      ss << "sv = \n" << svd.singularValues().format(python_fmt) << "\n";

      for (auto j : s) {
        ss << "j = " << j << "\n";
        ss << indent_block(
            1,
            string_format("points = %s\nweights = %s\n",
                          format_as_list(grid.cells[j].qr.weights).c_str(),
                          format_as_list(grid.cells[j].qr.points).c_str()));
      }

      LOG_WARN(ss.str());
    }

    return svd.rank() == A.cols();
  };

  auto s = conservative_stencil(grid, i_center, k, n_points);
  if (s.size() == n_points && is_good(s)) {
    return s;
  }

  s = less_conservative_stencil(grid, i_center, k, n_points);
  if (s.size() == n_points && is_good(s)) {
    return s;
  }

  return tryhard_stencil(grid, is_good, i_center, k, n_points);
}

std::vector<int_t>
central_stencil(const Grid &grid, int_t i_center, int_t n_points) {
  return region_based_stencil(grid, i_center, n_points, FullSphere());
}

} // namespace zisa
