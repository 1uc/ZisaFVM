#include <benchmark/benchmark.h>

#include <zisa/math/quadrature.hpp>
#include <zisa/reconstruction/cweno_ao.hpp>
#include <zisa/reconstruction/global_reconstruction.hpp>
#include <zisa/scenarios/euler_high_order.hpp>

namespace zisa {
namespace bm {
using namespace zisa::scenarios::euler_high_order;

static void zisa_global_reconstruction(benchmark::State &state) {
  auto grid = load_grid();
  auto global_reconstruction = make_global_reconstruction(grid);
  auto current_state = make_valid_initial_conditions(*grid);

  for (auto _ : state) {
    global_reconstruction->compute(*current_state);
  }
}

} // namespace bm
} // namespace zisa

static void bm_global_reconstruction(benchmark::State &state) {
  zisa::bm::zisa_global_reconstruction(state);
}

BENCHMARK(bm_global_reconstruction)->Unit(benchmark::kMicrosecond);
