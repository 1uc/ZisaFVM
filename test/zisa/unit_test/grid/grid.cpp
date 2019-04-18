#include <random>

#include <zisa/grid/grid.hpp>
#include <zisa/grid/grid_impl.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/barycentric.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/cartesian.hpp>
#include <zisa/math/poly2d.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/utils/to_string.hpp>

TEST_CASE("Grid; small_example", "[grid]") {

  zisa::int_t n_cells = 2;
  zisa::int_t n_vertices = 4;

  auto vertices = zisa::array<zisa::XYZ, 1>(zisa::shape_t<1>{n_vertices});
  vertices(0) = {0.0, 0.0, 0.0};
  vertices(1) = {1.0, 0.0, 0.0};
  vertices(2) = {1.0, 1.0, 0.0};
  vertices(3) = {0.0, 1.0, 0.0};

  auto vertex_indices
      = zisa::array<zisa::int_t, 2>(zisa::shape_t<2>{n_cells, 3ul});
  vertex_indices(0, 0) = 0;
  vertex_indices(0, 1) = 1;
  vertex_indices(0, 2) = 3;

  vertex_indices(1, 0) = 1;
  vertex_indices(1, 1) = 2;
  vertex_indices(1, 2) = 3;

  SECTION("compute_edges") {
    auto edges = zisa::compute_edges(vertex_indices);

    REQUIRE(edges[0].find(0) == edges[0].end());
    REQUIRE(edges[0].find(1) != edges[0].end());
    REQUIRE(edges[0].find(2) == edges[0].end());
    REQUIRE(edges[0].find(3) == edges[0].end());

    REQUIRE(edges[1].find(0) == edges[1].end());
    REQUIRE(edges[1].find(1) == edges[1].end());
    REQUIRE(edges[1].find(2) != edges[1].end());
    REQUIRE(edges[1].find(3) != edges[1].end());

    REQUIRE(edges[2].find(0) == edges[2].end());
    REQUIRE(edges[2].find(1) == edges[2].end());
    REQUIRE(edges[2].find(2) == edges[2].end());
    REQUIRE(edges[2].find(3) != edges[2].end());

    REQUIRE(edges[3].find(0) != edges[3].end());
    REQUIRE(edges[3].find(1) != edges[3].end());
    REQUIRE(edges[3].find(2) == edges[3].end());
    REQUIRE(edges[3].find(3) == edges[3].end());
  }

  SECTION("compute_neighbours") {
    auto neighbours = zisa::compute_neighbours(vertex_indices);

    REQUIRE(neighbours(0, 0) == zisa::magic_index_value);
    REQUIRE(neighbours(0, 1) == 1);
    REQUIRE(neighbours(0, 2) == zisa::magic_index_value);

    REQUIRE(neighbours(1, 0) == zisa::magic_index_value);
    REQUIRE(neighbours(1, 1) == zisa::magic_index_value);
    REQUIRE(neighbours(1, 2) == 0);
  }

  SECTION("compute_valid_neighbours") {
    auto neighbours = zisa::compute_neighbours(vertex_indices);
    auto is_valid = zisa::compute_valid_neighbours(neighbours);

    REQUIRE(!is_valid(0, 0));
    REQUIRE(is_valid(0, 1));
    REQUIRE(!is_valid(0, 2));

    REQUIRE(!is_valid(1, 0));
    REQUIRE(!is_valid(1, 1));
    REQUIRE(is_valid(1, 2));
  }

  SECTION("compute_edge_indices") {
    auto neighbours = zisa::compute_neighbours(vertex_indices);
    auto is_valid = zisa::compute_valid_neighbours(neighbours);
    auto edge_indices = zisa::compute_edge_indices(neighbours, is_valid);

    // Only interior (diagonal) edge.
    REQUIRE(edge_indices(0, 1) == 0);
    REQUIRE(edge_indices(1, 2) == 0);

    // Exterior edges.
    REQUIRE(edge_indices(0, 0) == 1);
    REQUIRE(edge_indices(0, 2) == 2);
    REQUIRE(edge_indices(1, 0) == 3);
    REQUIRE(edge_indices(1, 1) == 4);
  }

  SECTION("compute_normals") {
    auto neighbours = zisa::compute_neighbours(vertex_indices);
    auto is_valid = zisa::compute_valid_neighbours(neighbours);
    auto edge_indices = zisa::compute_edge_indices(neighbours, is_valid);
    auto normals = zisa::compute_normals(
        vertices, vertex_indices, neighbours, is_valid, edge_indices);

    REQUIRE(normals.shape(0) == 5);

    auto exact = zisa::XYZ{0.5 / std::sqrt(0.5), 0.5 / std::sqrt(0.5), 0.0};
    REQUIRE(zisa::almost_equal(normals(0), exact, 1e-12));
  }

  SECTION("compute_volumes") {
    auto volumes = compute_volumes(vertices, vertex_indices);

    REQUIRE(volumes.shape().size() == 1);
    REQUIRE(volumes.shape(0) == n_cells);
    REQUIRE(zisa::almost_equal(volumes(0), 0.5, 1e-12));
  }
}

TEST_CASE("Grid; sizes", "[grid]") {
  auto grid = zisa::load_gmsh("grids/dbg.msh");

  zisa::int_t n_cells = grid->n_cells;
  zisa::int_t n_vertices = grid->n_vertices;
  zisa::int_t n_edges = grid->n_edges;

  REQUIRE(grid->vertex_indices.shape().size() == 2);
  REQUIRE(grid->vertex_indices.shape(0) == n_cells);
  REQUIRE(grid->vertex_indices.shape(1) == 3);

  REQUIRE(grid->edge_indices.shape().size() == 2);
  REQUIRE(grid->edge_indices.shape(0) == n_cells);
  REQUIRE(grid->edge_indices.shape(1) == 3);

  REQUIRE(grid->neighbours.shape().size() == 2);
  REQUIRE(grid->neighbours.shape(0) == n_cells);
  REQUIRE(grid->neighbours.shape(1) == 3);

  REQUIRE(grid->is_valid.shape().size() == 2);
  REQUIRE(grid->is_valid.shape(0) == n_cells);
  REQUIRE(grid->is_valid.shape(1) == 3);

  REQUIRE(grid->vertices.shape().size() == 1);
  REQUIRE(grid->vertices.shape(0) == n_vertices);

  REQUIRE(grid->cell_centers.shape().size() == 1);
  REQUIRE(grid->cell_centers.shape(0) == n_cells);

  REQUIRE(grid->volumes.shape().size() == 1);
  REQUIRE(grid->volumes.shape(0) == n_cells);

  REQUIRE(grid->normals.shape().size() == 1);
  REQUIRE(grid->normals.shape(0) == n_edges);

  REQUIRE(grid->tangentials.shape().size() == 1);
  REQUIRE(grid->tangentials.shape(0) == n_edges);
}

TEST_CASE("Grid; moments", "[grid]") {
  auto grid = zisa::load_gmsh("grids/dbg.msh");

  auto check_moment
      = [](const zisa::Triangle &tri, int k, int l, double exact) {
          auto m = zisa::avg_moment(tri, k, l, 3);

          INFO(string_format(
              "[%d, %d] %e  !=  %e (%e) \n", k, l, m, exact, m - exact));
          REQUIRE(zisa::almost_equal(m, exact, 1e-8));
        };

  SECTION("avg_moment") {

    for (const auto &[i, tri] : zisa::triangles(*grid)) {

      check_moment(tri, 0, 0, 1.0);

      check_moment(tri, 1, 0, 0.0);
      check_moment(tri, 0, 1, 0.0);

      REQUIRE(zisa::abs(zisa::avg_moment(tri, 1, 2, 3)) < 1.0);
    }
  }

  SECTION("normalized_moments") {

    auto check_moment = [](const zisa::array<double, 1> &m, int k, int l) {
      auto ikl = zisa::poly_index(k, l);

      INFO(string_format("[%d, %d] %e  \n", k, l, m(ikl)));
      REQUIRE(zisa::abs(m(ikl)) < 2.0);
    };

    for (const auto &[i, tri] : zisa::triangles(*grid)) {
      auto m = zisa::normalized_moments(tri, 4, 3);

      for (int d = 0; d <= 3; ++d) {
        for (int k = 0; k <= d; ++k) {
          check_moment(m, k, d - k);
        }
      }
    }
  }

  SECTION("orientation normals") {
    for (const auto &[e, edge] : zisa::interior_edges(*grid)) {
      auto [iL, iR] = grid->left_right(e);

      zisa::XYZ xL = zisa::barycenter(grid->triangle(iL));
      zisa::XYZ xR = zisa::barycenter(grid->triangle(iR));

      auto n = edge.normal();

      // clang-format off
      INFO(string_format("[%d, (%d, %d)] xL = (%.3e, %.3e) xR = (%.3e, %.3e), "
                         "n = (%.3e, %.3e)",
                         e, iL, iR,
                         xL[0], xL[1],
                         xR[0], xR[1],
                         n[0], n[1]));
      // clang-format on
      REQUIRE(zisa::dot(n, xR - xL) > 0.0);
    }
  }
}

TEST_CASE("Grid; volume", "[grid]") {
  auto grid = zisa::load_gmsh("grids/convergence/unit_square_2.msh");

  REQUIRE(zisa::almost_equal(volume(*grid), 1.0, 1e-12));

  for (const auto &[i, tri] : triangles(*grid)) {
    REQUIRE(volume(tri) > 0.0);
  }
}

TEST_CASE("Grid; iterators", "[grid]") {
  auto grid = zisa::load_gmsh("grids/convergence/unit_square_2.msh");
  auto n_interior_edges = grid->n_interior_edges;
  auto n_cells = grid->n_cells;

  SECTION("interior_edges") {
    zisa::int_t count = 0;
    for (const auto &[e, edge] : interior_edges(*grid)) {
      ++count;
    }

    REQUIRE(count == n_interior_edges);
  }

  SECTION("triangles") {
    zisa::int_t count = 0;
    for (const auto &[i, tri] : triangles(*grid)) {
      ++count;
    }

    REQUIRE(count == n_cells);
  }
}

TEST_CASE("Grid; incidence", "[grid]") {
  auto grid = zisa::load_gmsh("grids/convergence/unit_square_1.msh");
  zisa::int_t n_cells = grid->n_cells;
  zisa::int_t max_neighbours = grid->max_neighbours;

  SECTION("(iL, iR) are neighbours") {
    for (const auto &[e, edge] : interior_edges(*grid)) {
      auto [iL, iR] = grid->left_right(e);

      REQUIRE((grid->neighbours(iL, 0) == iR || grid->neighbours(iL, 1) == iR
               || grid->neighbours(iL, 2) == iR));
    }
  }

  SECTION("edges of interior cells") {
    auto count = zisa::array<zisa::int_t, 1>(zisa::shape_t<1>{n_cells});

    for (const auto &[e, edge] : interior_edges(*grid)) {
      auto [iL, iR] = grid->left_right(e);

      count(iL) += 1;
      count(iR) += 1;
    }

    for (const auto &[i, tri] : triangles(*grid)) {
      zisa::int_t expected = 0;
      for (zisa::int_t k = 0; k < max_neighbours; ++k) {
        if (grid->is_valid(i, k)) {
          expected += 1;
        }
      }

      REQUIRE(expected == count(i));
    }
  }
}

TEST_CASE("Grid; serialize", "[grid]") {
  auto grid = zisa::load_gmsh("grids/small.msh");
  auto filename = std::string("__unit_tests--grid_to_hdf5.h5");

  auto writer = zisa::HDF5SerialWriter(filename);
  zisa::save(writer, *grid);
}

TEST_CASE("Grid; locate", "[grid]") {
  auto grid_names = std::vector<std::string>{};

  for (int i = 0; i < 3; ++i) {
    grid_names.push_back(
        string_format("grids/convergence/unit_square_%d.msh", i));
  }
  std::random_device rd;
  std::mt19937 gen(rd());

  for (const auto &grid_name : grid_names) {
    auto grid = zisa::load_gmsh(grid_name);

    auto n_cells = grid->n_cells;
    auto max_iter = grid->n_cells;

    std::uniform_int_distribution<zisa::int_t> i_guess_dis(0, n_cells - 1);
    std::uniform_real_distribution<double> x_dis(0.05, 0.95);
    for (auto [i, tri] : triangles(*grid)) {
      auto i_guess = i_guess_dis(gen);

      auto eta = x_dis(gen);
      auto zeta = x_dis(gen) * (1.0 - eta);
      auto x = zisa::coord(tri, zisa::Barycentric{eta, zeta, 1.0 - eta - zeta});

      auto i_cell = zisa::locate(*grid, x, i_guess, max_iter);

      INFO(string_format("i_guess = %d, tri = %s, x = %s",
                         i_guess,
                         zisa::to_string(tri).c_str(),
                         zisa::to_string(x).c_str()));
      REQUIRE(*i_cell == i);
    }
  }
}
