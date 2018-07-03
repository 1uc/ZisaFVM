#include <catch/catch.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/grid/grid_impl.hpp>
#include <zisa/math/basic_functions.hpp>
#include <zisa/math/cartesian.hpp>

TEST_CASE("Grid; small_example", "[grid]") {
  auto vertices = zisa::array<zisa::XY, 1>(zisa::shape_t<1>{4ul});
  vertices(0) = {0.0, 0.0};
  vertices(1) = {1.0, 0.0};
  vertices(2) = {1.0, 1.0};
  vertices(3) = {0.0, 1.0};

  auto vertex_indices = zisa::array<zisa::int_t, 2>(zisa::shape_t<2>{2ul, 3ul});
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

    REQUIRE(edge_indices(0, 0) == zisa::magic_index_value);
    REQUIRE(edge_indices(0, 1) == 0);
    REQUIRE(edge_indices(0, 2) == zisa::magic_index_value);

    REQUIRE(edge_indices(1, 0) == zisa::magic_index_value);
    REQUIRE(edge_indices(1, 1) == zisa::magic_index_value);
    REQUIRE(edge_indices(1, 2) == 0);
  }

  SECTION("compute_normals") {
    auto neighbours = zisa::compute_neighbours(vertex_indices);
    auto is_valid = zisa::compute_valid_neighbours(neighbours);
    auto edge_indices = zisa::compute_edge_indices(neighbours, is_valid);
    auto normals = zisa::compute_normals(
        vertices, vertex_indices, neighbours, is_valid, edge_indices);

    REQUIRE(normals.shape(0) == 1);

    auto exact = zisa::XY{0.5 / std::sqrt(0.5), 0.5 / std::sqrt(0.5)};
    REQUIRE(zisa::almost_equal(normals(0), exact, 1e-12));
  }

  SECTION("compute_volumes") {
    auto volumes = compute_volumes(vertices, vertex_indices);
    REQUIRE(zisa::almost_equal(volumes(0), 0.5, 1e-12));
  }
}
