#include <benchmark/benchmark.h>

#include <random>
#include <zisa/math/poly2d.hpp>

static void bm_poly2d_eval(benchmark::State &state) {
  std::random_device rd;
  std::mt19937 e(rd());
  std::uniform_real_distribution<double> dist(0, 1.0);

  auto rand = [&e, &dist]() { return dist(e); };

  // clang-format off
  auto p = zisa::Poly2D<4, 1>(
                           {rand(),
                            rand(), rand(),
                            rand(), rand(), rand(),
                            rand(), rand(), rand(), rand(),
                            rand(), rand(), rand(), rand(), rand()},

                           {0.0,
                            0.0, 0.0,
                            rand(), rand(), rand(),
                            rand(), rand(), rand(), rand(),
                            rand(), rand(), rand(), rand(), rand()},

                           zisa::XYZ{rand(), rand(), rand()},
                           0.234);
  // clang-format on

  auto x = zisa::XYZ{rand(), rand(), 0.0};

  for (auto _ : state) {
    benchmark::DoNotOptimize(p(x));
  }
}

BENCHMARK(bm_poly2d_eval);
