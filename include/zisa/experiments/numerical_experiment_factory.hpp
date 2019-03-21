#ifndef NUMERICAL_EXPERIMENT_FACTORY_H_8T8A7
#define NUMERICAL_EXPERIMENT_FACTORY_H_8T8A7

#include <zisa/cli/input_parameters.hpp>
#include <zisa/config.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/experiments/numerical_experiment_factory.hpp>

namespace zisa {

class NumericalExperimentFactory {
private:
  using return_type = std::unique_ptr<zisa::NumericalExperiment>;
  using argument_type = InputParameters;
  using callback_type = std::function<return_type(const argument_type &)>;
  using key_type = std::string;

private:
  NumericalExperimentFactory() = default;
  NumericalExperimentFactory(const NumericalExperimentFactory &) = delete;
  NumericalExperimentFactory(NumericalExperimentFactory &&) = delete;

  ~NumericalExperimentFactory() = default;

public:
  static NumericalExperimentFactory &instance();
  void add(const key_type &key, const callback_type &callback);
  return_type make(const key_type &key, const argument_type &args) const;

protected:
  bool is_good_key(const key_type &key) const;

private:
  std::map<key_type, callback_type> registry_;
};

void register_experiment(
    const std::string &key,
    const std::function<
        std::unique_ptr<NumericalExperiment>(const InputParameters &)> &f);

std::unique_ptr<zisa::NumericalExperiment>
make_experiment(const InputParameters &params);

} // namespace zisa

#endif
