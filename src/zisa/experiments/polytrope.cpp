#include <zisa/experiments/ic/polytrope_ic.hpp>
#include <zisa/experiments/polytrope.hpp>

namespace zisa {

std::shared_ptr<AllVariables> Polytrope::choose_initial_conditions() {
  auto dims = choose_all_variable_dims();
  auto all_variables = std::make_shared<AllVariables>(dims);

  auto qr = choose_volume_rule();
  const auto &eos = euler.eos;

  auto ic_ = PolytropeIC(euler);
  auto ic = [&eos, &ic_](const auto &x) { return eos.cvars(ic_(x)); };

  auto &u0 = all_variables->cvars;
  for (auto &&[i, tri] : triangles(*grid)) {
    u0(i) = average(qr, ic, tri);
  }

  return all_variables;
}

} // namespace zisa
