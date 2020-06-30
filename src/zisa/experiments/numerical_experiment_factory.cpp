#include <zisa/experiments/numerical_experiment_factory.hpp>

// include experiment headers.
#include <zisa/experiments/polytrope.hpp>
#include <zisa/experiments/rayleigh_taylor.hpp>
#include <zisa/experiments/smooth_bubble.hpp>

#ifdef ZISA_HAS_MPI
#include <zisa/experiments/mpi_numerical_experiment.hpp>
#endif

namespace zisa {

class NumericalExperimentFactory {
private:
  using return_type = std::unique_ptr<zisa::NumericalExperiment>;
  using argument_type = InputParameters;
  using callback_type = std::function<return_type(const argument_type &)>;
  using key_type = std::string;

public:
  return_type make(const key_type &key, const argument_type &args) const;

  template <class Experiment>
  void register_simple(const key_type &key);
  void register_generic(const key_type &key, const callback_type &callback);

protected:
  bool is_good_key(const key_type &key) const;

private:
  std::map<key_type, callback_type> registry_;
};

void NumericalExperimentFactory::register_generic(
    const key_type &key, const callback_type &callback) {

  assert(is_good_key(key));
  registry_[key] = callback;
}

template <class Experiment>
void NumericalExperimentFactory::register_simple(const key_type &key) {
#if ZISA_HAS_MPI == 1
  register_generic(
      key,
      [](const InputParameters &params)
          -> std::unique_ptr<NumericalExperiment> {
        if (is_mpi(params)) {
          return std::make_unique<MPINumericalExperiment<Experiment>>(params);
        } else {
          return std::make_unique<Experiment>(params);
        }
      });
#else
  register_generic(key, [](const InputParameters &params) {
    return std::make_unique<Experiment>(params);
  });
#endif
}

bool NumericalExperimentFactory::is_good_key(const key_type &key) const {
  return !key.empty() && registry_.find(key) == registry_.cend();
}

auto NumericalExperimentFactory::make(const key_type &key,
                                      const argument_type &args) const
    -> return_type {

  LOG_ERR_IF(registry_.find(key) == registry_.cend(),
             string_format("Unknown numerical experiment. [%s]", key.c_str()));

  return registry_.at(key)(args);
}

static NumericalExperimentFactory make_factory() {
  NumericalExperimentFactory factory;

  factory.register_simple<SmoothBubble>("smooth_bubble");
  factory.register_simple<Polytrope>("scaling_experiment");
  factory.register_simple<Polytrope>("gaussian_bump");
  factory.register_simple<Polytrope>("gaussian_bump_3d");
  factory.register_simple<RayleighTaylor>("rayleigh_taylor");

  return factory;
}

std::unique_ptr<NumericalExperiment>
make_experiment(const InputParameters &params) {
  static auto factory = make_factory();

  std::string exp = params["experiment"]["name"];
  return factory.make(exp, params);
}

} // namespace zisa
