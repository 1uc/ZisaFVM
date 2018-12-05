
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/experiments/shock_bubble.hpp>

namespace zisa {

std::unique_ptr<NumericalExperiment>
make_experiment(const InputParameters &params) {

  std::string exp = params["experiment"];
  if (exp == "shock_bubble") {
    return std::make_unique<ShockBubble>(params);
  }

  LOG_ERR(string_format("Unknown numerical experiment. [%s]", exp.c_str()));
}

template <>
ConstantGravityRadial
make_gravity<ConstantGravityRadial>(const InputParameters &params) {
  LOG_ERR_IF(!has_key(params, "euler"), "Missing section 'euler'");
  LOG_ERR_IF(!has_key(params["euler"], "gravity"), "Missing section 'euler/gravity'");
  LOG_ERR_IF(!has_key(params["euler"]["gravity"], "g"), "Missing entry 'euler/gravity/g'");

  return ConstantGravityRadial(double(params["euler"]["gravity"]["g"]));
}

} // namespace zisa
