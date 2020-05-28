#include <zisa/grid/grid.hpp>
#include <zisa/memory/tree.hpp>
#include <zisa/testing/testing_framework.hpp>
#include <zisa/unit_test/grid/test_grid_factory.hpp>

TEST_CASE("Tree; basic API", "[.][memory][tree]") {
  auto grid = zisa::load_grid(zisa::TestGridFactory::small());

  auto tree = zisa::Tree<int, 4>{};
  auto f = [&grid](zisa::int_t i, zisa::int_t j) {
    return zisa::norm(grid->cell_centers(i) - grid->cell_centers(j));
  };

  REQUIRE(tree.empty());

  for (zisa::int_t i = 0; i < 4; ++i) {
    tree.insert(i, f);
  }

  REQUIRE(!tree.empty());

  for (zisa::int_t i = 0; i < 4; ++i) {
    auto look_up = zisa::int_t(i);
    auto found
        = tree.locate([&f, &look_up](zisa::int_t i) { return f(look_up, i); });

    REQUIRE(look_up == found);
  }
}
