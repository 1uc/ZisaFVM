#ifndef NUMERICAL_EXPERIMENT_FACTORY_H_8T8A7
#define NUMERICAL_EXPERIMENT_FACTORY_H_8T8A7

#include <zisa/cli/input_parameters.hpp>
#include <zisa/config.hpp>
#include <zisa/experiments/numerical_experiment.hpp>
#include <zisa/experiments/numerical_experiment_factory.hpp>

namespace zisa {

std::unique_ptr<zisa::NumericalExperiment>
make_experiment(const InputParameters &params);

} // namespace zisa

#endif
