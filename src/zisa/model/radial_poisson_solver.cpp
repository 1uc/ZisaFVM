
#include <zisa/math/spherical_shell.hpp>
#include <zisa/model/euler_variables.hpp>
#include <zisa/model/radial_poisson_solver.hpp>

namespace zisa {

RadialPoissonSolver::RadialPoissonSolver(std::shared_ptr<Grid> grid,
                                         array<array<int_t, 1>, 1> cell_indices,
                                         double gravitational_constant)
    : grid(std::move(grid)),
      cell_indices(std::move(cell_indices)),
      G(gravitational_constant) {}

void RadialPoissonSolver::update(RadialGravity &gravity,
                                 const AllVariables &all_variables) const {

  update_gravity(
      gravity, [&all_variables](int_t i) { return all_variables.cvars(i, 0); });
}

template <class Rho>
void RadialPoissonSolver::update_gravity(RadialGravity &gravity,
                                         const Rho &rho) const {

  auto &radii = gravity.radius_array();
  auto &phi = gravity.phi_array();

  double rho_center = rho(cell_indices[0][0]);
  double m_enc = 0.0;
  phi[0] = 0.0;
  phi[1] = 2.0 / 3.0 * G * zisa::pi * rho_center * zisa::pow<2>(radii[1]);

  int_t n_radii = radii.size();
  for (int_t l = 1; l < n_radii - 1; ++l) {
    m_enc += layer_mass(rho, radii, l - 1);

    double dr = radii[l + 1] - radii[l];
    LOG_ERR_IF(dr <= 0.0, "non-positive dr");

    double r2 = zisa::pow<2>(0.5 * (radii[l + 1] + radii[l]));
    phi[l + 1] = phi[l] + dr * G * m_enc / r2;

    LOG_ERR_IF(phi[l + 1] - phi[l] <= 0.0, "non-positve increase.");
  }

  phi[n_radii - 1] = phi[n_radii - 2] * (1.0 + 1e-12);
}

template <class Rho>
double RadialPoissonSolver::layer_mass(const Rho &rho,
                                       const array<double, 1> &radii,
                                       int_t layer) const {
  int_t n_cells_per_layer = cell_indices[layer].size();

  if (n_cells_per_layer == 0) {
    //    LOG_WARN(string_format("In layer %d there are zero cells.", layer));
    return 0.0;
  }

  double rho_avg = 0.0;
  double volume_shell = 0.0;

#if ZISA_HAS_OPENMP == 1
#pragma omp parallel for reduction(+ : rho_avg)
#endif
  for (int_t i_loc = 0; i_loc < n_cells_per_layer; ++i_loc) {
    int_t i = cell_indices[layer][i_loc];

    double vol = grid->volumes(i);
    rho_avg += vol * rho(i);
    volume_shell += vol;
  }
  rho_avg = rho_avg / volume_shell;

  double r_lower = zisa::avg(radii[layer], radii[layer + 1]);
  double r_outer = zisa::avg(radii[layer + 1], radii[layer + 2]);

  if (layer == 0) {
    r_lower = radii[0];
  }

  if (layer == radii.size() - 3) {
    r_outer = radii[radii.size() - 1];
  }

  return rho_avg * volume(SphericalShell{r_lower, r_outer});
}

array<double, 1> make_radial_bins(const Grid &grid,
                                  double r_inner,
                                  double r_outer,
                                  double rel_layer_width) {
  std::vector<double> radii_;
  radii_.reserve(int_t(zisa::sqrt(double(grid.n_cells))));
  radii_.push_back(r_inner);

  auto locate_next_cell = [&grid, &radii_]() {
    return zisa::locate(grid, XYZ(radii_.back() * XYZ::unit_vector(0)));
  };

  std::optional<int_t> i_cell;
  while ((i_cell = locate_next_cell()) != std::nullopt) {
    double dr = 2.0 * circum_radius(grid.triangle(i_cell.value()));
    radii_.push_back(radii_.back() + rel_layer_width * dr);
  }
  radii_.back() = r_outer;

  auto radii = zisa::array<double, 1>(radii_.size());
  std::copy(radii_.cbegin(), radii_.cend(), radii.begin());

  return radii;
}

void save(HierarchicalWriter &writer, const RadialPoissonSolver &solver) {
  writer.open_group("poisson_solver");
  writer.write_string("RadialPoissonSolver", "name");
  writer.write_scalar(solver.G, "gravitational_constant");

  writer.open_group("cell_indices");

  const auto &cell_indices = solver.cell_indices;
  auto n_layers = cell_indices.size();
  writer.write_scalar(n_layers, "size");

  for (int_t i = 0; i < n_layers; ++i) {
    auto tag = string_format("%d", i);
    save(writer, cell_indices[i], tag);
  }
  writer.close_group();
  writer.close_group();
}

RadialPoissonSolver
RadialPoissonSolver::load(HierarchicalReader &reader,
                          const std::shared_ptr<Grid> &grid) {

  reader.open_group("poisson_solver");

  auto name = reader.read_string("name");
  auto G = reader.read_scalar<double>("gravitational_constant");

  reader.open_group("cell_indices");

  auto n_layers = reader.read_scalar<int_t>("size");
  auto cell_indices = array<array<int_t, 1>, 1>(shape_t<1>{n_layers});

  for (int_t i = 0; i < n_layers; ++i) {
    auto tag = string_format("%d", i);
    cell_indices[i] = array<int_t, 1>::load(reader, tag);
  }
  reader.close_group();
  reader.close_group();

  return RadialPoissonSolver(grid, std::move(cell_indices), G);
}

array<double, 1>
make_radial_bins(const Grid &grid, double r_outer, double rel_layer_width) {
  return make_radial_bins(grid, 0.0, r_outer, rel_layer_width);
}

class HalfRadii {
public:
  class Iterator {
  public:
    explicit Iterator(const array<double, 1> &radii, int_t i)
        : i(i), radii_(radii) {}

    void operator++() { ++i; }

    double operator*() const {
      auto n_radii = radii_.size();
      if (i == 0) {
        return radii_[0];
      } else if (i < n_radii - 1) {
        return 0.5 * (radii_[i] + radii_[i + 1]);
      } else {
        return radii_[n_radii - 1];
      }
    }

    bool operator==(const Iterator &other) const { return i == other.i; }

    bool operator!=(const Iterator &other) const { return !((*this) == other); }

    int_t index() const { return i; }

  private:
    int_t i;
    const array<double, 1> &radii_;
  };

  explicit HalfRadii(const array<double, 1> &radii) : radii(radii) {}

  Iterator begin() const { return Iterator(radii, 0); }
  Iterator end() const { return Iterator(radii, radii.size() - 1); }

private:
  const array<double, 1> &radii;
};
}

namespace std {
template <>
struct iterator_traits<zisa::HalfRadii::Iterator> {
  using iterator_category = std::forward_iterator_tag;
};
}

namespace zisa {
array<array<int_t, 1>, 1>
make_cell_indices_bins(const Grid &grid, const array<double, 1> &radii) {

  auto half_radii = HalfRadii(radii);

  int_t n_layers = radii.size() - 2;
  auto cell_indices_ = std::vector<std::vector<int_t>>(n_layers);
  for (auto &ci : cell_indices_) {
    ci.reserve(100);
  }

  auto max_neighbours = grid.max_neighbours;

  for (const auto &[i, tri] : triangles(grid)) {
    //    double r = zisa::norm(barycenter(tri));
    //    auto iter = std::find_if(half_radii.begin(),
    //                             half_radii.end(),
    //                             [r](double r_test) { return r < r_test; });
    //    auto l = iter.index() - 1;
    //    cell_indices_[l].push_back(i);

    for (int_t k = 0; k < max_neighbours; ++k) {

      double r = zisa::norm(grid.vertex(i, k));
      auto iter = std::find_if(half_radii.begin(),
                               half_radii.end(),
                               [r](double r_test) { return r < r_test; });

      auto l = iter.index() - 1;
      l = zisa::min(l, cell_indices_.size() - 1);
      cell_indices_[l].push_back(i);
    }
  }

  auto cell_indices = array<array<int_t, 1>, 1>(n_layers);
  for (int_t l = 0; l < n_layers; ++l) {
    LOG_ERR_IF(cell_indices_[l].size() == 0, string_format("fail. [%d]", l));

    cell_indices[l] = array<int_t, 1>(cell_indices_[l].size());
    std::copy(cell_indices_[l].cbegin(),
              cell_indices_[l].cend(),
              cell_indices[l].begin());
  }

  return cell_indices;
}

std::pair<RadialGravity, std::shared_ptr<RadialPoissonSolver>>
make_radial_poisson_solver(const std::shared_ptr<Grid> &grid,
                           double gravitational_constant) {

  double r_outer = 0.0;
  for (int_t i = 0; i < grid->n_vertices; ++i) {
    r_outer = zisa::max(r_outer, zisa::norm(grid->vertices[i]));
  }

  auto radii = make_radial_bins(*grid, r_outer, 1.0);
  auto phi = array<double, 1>(radii.shape());
  auto gravity = RadialGravity(std::move(radii), std::move(phi));

  auto ps = make_radial_poisson_solver(
      gravity.radius_array(), grid, gravitational_constant);
  return {gravity, ps};
}

std::shared_ptr<RadialPoissonSolver>
make_radial_poisson_solver(const array<double, 1> &radii,
                           const std::shared_ptr<Grid> &grid,
                           double gravitational_constant) {

  auto cell_indices = make_cell_indices_bins(*grid, radii);

  return std::make_shared<RadialPoissonSolver>(
      grid, std::move(cell_indices), gravitational_constant);
}

// double CrudeRadialPoissonSolver::layer_mass(int_t layer) const {
//  const auto &ci = cell_indices[layer];
//  auto &radii = euler->gravity.radius_array();
//
//  int_t n_quad_points = ci.size();
//  double rho_bar = 0.0;
//  for (int_t i = 0; i < n_quad_points; ++i) {
//    auto x = XYZ(zisa::avg(radii[layer], radii[layer + 1]) * quad_points[i]);
//    rho_bar += (*grc)(ci[i], x)[0];
//  }
//
//  rho_bar /= double(n_quad_points);
//
//  return rho_bar * volume(SphericalShell{radii[layer], radii[layer + 1]});
//}
//
// void CrudeRadialPoissonSolver::compute(const SimulationClock &,
//                                       AllVariables &) {
//  auto &phi = euler->gravity.phi_array();
//  auto &radii = euler->gravity.radius_array();
//
//  int_t n_layers = radii.size() - 1;
//
//  double m_enc = 0.0;
//  for (int_t l = 0; l < n_layers; ++l) {
//    m_enc += layer_mass(l);
//    double dr = radii[l + 1] - radii[l];
//    double r2 = zisa::pow<2>(radii[l + 1]);
//    phi[l + 1] = phi[l] + dr * G * m_enc / r2;
//  }
//}
//
// std::shared_ptr<CrudeRadialPoissonSolver> make_crude_radial_poisson_solver(
//    const std::shared_ptr<Euler<JankaEOS, RadialGravity>> &euler,
//    const std::shared_ptr<GlobalReconstruction<euler_var_t>> &grc,
//    const Grid &grid,
//    double gravitational_constant,
//    int_t n_quad_points) {
//
//  const auto &radii = euler->gravity.radius_array();
//
//  auto quad_points = array<XYZ, 1>(n_quad_points);
//  for (int_t i = 0; i < n_quad_points; ++i) {
//    quad_points[i]
//        = XYZ::r_hat(double(i) * 2.0 * zisa::pi / double(n_quad_points));
//  }
//
//  int_t n_layers = radii.size() - 1;
//  auto cell_indices = array<array<int_t, 1>, 1>(n_layers);
//
//  int_t max_iter = grid.n_cells;
//  int_t i_guess = zisa::locate(grid, XYZ::zeros()).value();
//
//  for (int_t l = 0; l < n_layers; ++l) {
//    cell_indices[l] = array<int_t, 1>(n_quad_points);
//
//    for (int_t i = 0; i < n_quad_points; ++i) {
//      auto x = XYZ(zisa::avg(radii[l], radii[l + 1]) * quad_points[i]);
//      std::optional<int_t> o = zisa::locate(grid, x, i_guess, max_iter);
//      cell_indices[l][i] = o.value();
//      i_guess = o.value();
//    }
//  }
//
//  return std::make_shared<CrudeRadialPoissonSolver>(euler,
//                                                    grc,
//                                                    std::move(cell_indices),
//                                                    std::move(quad_points),
//                                                    gravitational_constant);
//}

}
