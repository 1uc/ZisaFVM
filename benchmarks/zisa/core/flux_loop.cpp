#include <benchmark/benchmark.h>

#include <zisa/scenarios/euler_high_order.hpp>

namespace zisa {
namespace bm {
using namespace zisa::scenarios::euler_high_order;

static void zisa_flux_loop(benchmark::State &state) {
  auto grid = load_grid();
  auto roc = make_flux_loop(grid);

  auto tendency = make_all_variables(grid->n_cells);
  auto current_state = make_valid_initial_conditions(*grid);

  for (auto _ : state) {
    roc->compute(*tendency, *current_state, /* t = */ 0.0);
  }
}

static void zisa_source_loop(benchmark::State &state) {
  auto grid = load_grid();
  auto roc = make_source_loop(grid);

  auto tendency = make_all_variables(grid->n_cells);
  auto current_state = make_valid_initial_conditions(*grid);

  for (auto _ : state) {
    roc->compute(*tendency, *current_state, /* t = */ 0.0);
  }
}

static void zisa_zero_rate_of_change(benchmark::State &state) {
  auto grid = load_grid();
  auto roc = make_zero_rate_of_change();

  auto tendency = make_all_variables(grid->n_cells);
  auto current_state = make_valid_initial_conditions(*grid);

  for (auto _ : state) {
    roc->compute(*tendency, *current_state, /* t = */ 0.0);
  }
}

static void zisa_flux_bc(benchmark::State &state) {
  auto grid = load_grid();
  auto roc = make_flux_bc(grid);

  auto tendency = make_all_variables(grid->n_cells);
  auto current_state = make_valid_initial_conditions(*grid);

  for (auto _ : state) {
    roc->compute(*tendency, *current_state, /* t = */ 0.0);
  }
}

} // namespace bm
} // namespace zisa

static void bm_flux_loop(benchmark::State &state) {
  zisa::bm::zisa_flux_loop(state);
}

static void bm_source_loop(benchmark::State &state) {
  zisa::bm::zisa_source_loop(state);
}

static void bm_flux_bc(benchmark::State &state) {
  zisa::bm::zisa_flux_bc(state);
}

static void bm_zero_rate_of_change(benchmark::State &state) {
  zisa::bm::zisa_zero_rate_of_change(state);
}

BENCHMARK(bm_flux_loop)->Unit(benchmark::kMicrosecond);
BENCHMARK(bm_source_loop)->Unit(benchmark::kMicrosecond);
BENCHMARK(bm_flux_bc)->Unit(benchmark::kNanosecond);
BENCHMARK(bm_zero_rate_of_change)->Unit(benchmark::kMicrosecond);
