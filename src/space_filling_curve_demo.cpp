#include <zisa/grid/grid.hpp>
#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/math/space_filling_curve.hpp>

namespace zisa {

void space_filling_curve_demo() {

  auto grid = zisa::load_grid("grids/gaussian_bump-2.msh.h5");
  auto n_cells = grid->n_cells;

  std::vector<double> xs(n_cells);
  std::vector<double> ys(n_cells);

  for (int_t i = 0; i < n_cells; ++i) {
    xs[i] = grid->cell_centers[i][0];
    ys[i] = grid->cell_centers[i][1];
  }

  auto [x_min_, x_max_] = std::minmax_element(xs.begin(), xs.end());
  auto [y_min_, y_max_] = std::minmax_element(ys.begin(), ys.end());

  double x_min = *x_min_, x_max = *x_max_;
  double y_min = *y_min_, y_max = *y_max_;

  for (int_t i = 0; i < n_cells; ++i) {
    xs[i] = (xs[i] - x_min) / (x_max - x_min + 1e-10 * x_max);
    ys[i] = (ys[i] - y_min) / (y_max - y_min + 1e-10 * y_max);
  }

  std::vector<unsigned long long> sfc_indices(n_cells);
  std::vector<int_t> permutation(n_cells);

  for (int_t i = 0; i < n_cells; ++i) {
    sfc_indices[i] = hilbert_index<32>(xs[i], ys[i]).to_ullong();
    permutation[i] = i;
  }

  std::sort(
      permutation.begin(), permutation.end(), [&sfc_indices](int_t i, int_t j) {
        return sfc_indices[i] < sfc_indices[j];
      });

  auto n_parts = std::vector<int_t>{2, 4, 8, 16, 32, 64, 128, 256};

  {
    auto writer = HDF5SerialWriter("sfc_grid.h5");
    save(writer, *grid);
  }

  auto writer = HDF5SerialWriter("sfc_parts.h5");
  auto part = array<int_t, 1>(n_cells);
  for (auto np : n_parts) {
    auto block_size = n_cells / np;
    for (int_t i = 0; i < n_cells; ++i) {
      part[permutation[i]] = i / block_size;
    }

    save(writer, part, string_format("%04d", np));
  }
}

}

int main() {
  zisa::space_filling_curve_demo();
  return EXIT_SUCCESS;
}