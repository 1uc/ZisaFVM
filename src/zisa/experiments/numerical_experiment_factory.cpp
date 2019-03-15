
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/experiments/polytrope.hpp>
#include <zisa/experiments/shock_bubble.hpp>

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

  if (exp == "gaussian_bump") {
    return std::make_unique<Polytrope>(params);
  }

  LOG_ERR(string_format("Unknown numerical experiment. [%s]", exp.c_str()));
}

} // namespace zisa
