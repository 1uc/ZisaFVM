#include <memory>

#include <zisa/grid/grid.hpp>
#include <zisa/grid/point_locator.hpp>
#include <zisa/math/bounding_box.hpp>
#include <zisa/memory/tree.hpp>

namespace zisa {

PointLocator::PointLocator(std::shared_ptr<Grid> grid,
                           std::shared_ptr<Tree<int_t, 4>> &tree)
    : grid(std::move(grid)), tree(std::move(tree)) {}

int_t PointLocator::locate(const XYZ &x) const {
  auto i_guess = tree->locate(
          [this, &x](int_t i) { return zisa::norm(grid->cell_centers(i) - x); });

  auto i_cell = zisa::locate(*grid, x, i_guess, /* max_iter = */ 10);

  LOG_ERR_IF(!i_cell, "Could not find a cell containing the point.");

  return *i_cell;
}

std::shared_ptr<PointLocator>
make_point_locator(const std::shared_ptr<Grid> &grid) {
  auto tree = std::make_shared<Tree<int_t, 4>>();

  auto n_cells = grid->n_cells;
  auto bb = bounding_box(*grid);

  auto coord = [&bb](const XYZ rel) {
    const auto &[m, M] = bb;
    return XYZ(m + rel * (M - m));
  };

  auto f = [&grid = *grid](int_t i, int_t j) {
    return zisa::norm(grid.cell_centers(i) - grid.cell_centers(j));
  };

  int_t count = 0;
  for (int l = 1; count < n_cells; ++l) {
    auto n_points = 1 << l;
    for (int i = 0; i < n_points; ++i) {
      double dx = 1.0 / n_points;
      double x_rel = (i + 0.5) * dx;

      for (int j = 0; j < n_points; ++j) {
        double y_rel = (j + 0.5) * dx;

        auto x = coord(XYZ{x_rel, y_rel, 0});

        int_t i_guess = 0;
        if (!tree->empty()) {
          i_guess = tree->locate([&x, &grid = *grid](int_t i) {
            return zisa::norm(grid.cell_centers(i) - x);
          });
        }

        auto i_cell = locate(*grid, x, i_guess, grid->n_cells);
        ++count;

        if (i_cell) {
          tree->insert(*i_cell, f);
        }
      }
    }
  }

  return std::make_shared<PointLocator>(grid, tree);
}
}
