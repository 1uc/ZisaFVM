#ifndef ZISA_POINT_LOCATOR_HPP_234IH
#define ZISA_POINT_LOCATOR_HPP_234IH

#include <memory>

#include <zisa/config.hpp>
#include <zisa/grid/grid.hpp>
#include <zisa/memory/tree.hpp>

namespace zisa {

/// Locate a point in a *convex* domain.
class PointLocator {
public:
  explicit PointLocator(std::shared_ptr<Grid> grid,
                        std::shared_ptr<Tree<int_t, 4>> &tree);

  int_t locate(const XYZ &x) const {
    auto i_guess = tree->locate(
        [this, &x](int_t i) { return zisa::norm(grid->cell_centers(i) - x); });

    auto i_cell = zisa::locate(*grid, x, i_guess, /* max_iter = */ 10);

    LOG_ERR_IF(!i_cell, "Could not find a cell containing the point.");

    return *i_cell;
  }

private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<Tree<int_t, 4>> tree;
};

std::shared_ptr<PointLocator>
make_point_locator(const std::shared_ptr<Grid> &grid);

}
#endif // ZISA_POINT_LOCATOR_HPP
