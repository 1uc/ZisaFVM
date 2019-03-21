
#include <zisa/experiments/numerical_experiment_factory.hpp>
#include <zisa/experiments/polytrope.hpp>
#include <zisa/experiments/shock_bubble.hpp>

namespace zisa {

NumericalExperimentFactory &NumericalExperimentFactory::instance() {
  static NumericalExperimentFactory instance_;
  return instance_;
}

void NumericalExperimentFactory::add(const key_type &key,
                                     const callback_type &callback) {
  assert(is_good_key(key));
  registry_[key] = callback;
}

bool NumericalExperimentFactory::is_good_key(const key_type &key) const {
  return !key.empty() && registry_.find(key) == registry_.cend();
}

auto NumericalExperimentFactory::make(const key_type &key,
                                      const argument_type &args) const
    -> return_type {

  return registry_.at(key)(args);
}

std::unique_ptr<NumericalExperiment>
make_experiment(const InputParameters &params) {
  const auto &factory = NumericalExperimentFactory::instance();

  std::string exp = params["experiment"]["name"];
  return factory.make(exp, params);
}

void register_experiment(
    const std::string &key,
    const std::function<
        std::unique_ptr<NumericalExperiment>(const InputParameters &)> &f) {

  auto &instance = NumericalExperimentFactory::instance();
  instance.add(key, f);
}

} // namespace zisa
