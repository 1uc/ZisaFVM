#include <zisa/math/permutation.hpp>

namespace zisa {

int_t count_cycle_length(const array_const_view<int_t, 1> &sigma,
                         int_t i_start_cycle) {
  int_t i_next = sigma(i_start_cycle);
  int_t n_current = 1;
  while (i_next != i_start_cycle) {
    n_current += 1;
    i_next = sigma(i_next);
  }

  return n_current;
}

Permutation factor_permutation(const array_const_view<int_t, 1> &sigma) {
  std::vector<Cycle> cycles;
  cycles.reserve(100);

  auto is_processed = std::vector<bool>(sigma.size(), false);
  int_t i_next_cycle = 0;
  int_t n_cells = sigma.size();

  while (i_next_cycle < n_cells) {
    int_t i_start_cycle = i_next_cycle;
    int_t i_current = i_start_cycle;

    int_t n_current = count_cycle_length(sigma, i_start_cycle);
    cycles.emplace_back(n_current);
    auto &cycle = cycles.back();

    for (int_t i = 0; i < n_current; ++i) {
      cycle(i) = i_current;
      is_processed[i_current] = true;

      if (i_next_cycle == i_current) {
        while (i_next_cycle < n_cells && is_processed[i_next_cycle]) {
          ++i_next_cycle;
        }
      }

      i_current = sigma(i_current);
    }
  }

  int_t n_cycles = cycles.size();
  auto cycles_array = array<Cycle, 1>(n_cycles);
  for (int_t i = 0; i < n_cycles; ++i) {
    cycles_array(i) = std::move(cycles[i]);
  }

  return Permutation{std::move(cycles_array)};
}
}