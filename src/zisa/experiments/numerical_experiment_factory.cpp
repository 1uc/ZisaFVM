
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/experiments/shock_bubble.hpp>
#include <zisa/experiments/polytrope.hpp>

namespace zisa {

std::unique_ptr<NumericalExperiment>
make_experiment(const InputParameters &params) {

  std::string exp = params["experiment"]["name"];
  if (exp == "shock_bubble") {
    return std::make_unique<ShockBubble>(params);
  }

  if (exp == "polytrope") {
    return std::make_unique<Polytrope>(params);
  }


  LOG_ERR(string_format("Unknown numerical experiment. [%s]", exp.c_str()));
}

} // namespace zisa
