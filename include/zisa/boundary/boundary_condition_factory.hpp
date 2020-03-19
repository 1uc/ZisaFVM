#ifndef ZISA_BOUNDARY_CONDITION_FACTORY_HPP_CXLZKL
#define ZISA_BOUNDARY_CONDITION_FACTORY_HPP_CXLZKL

#include <zisa/config.hpp>

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/boundary/frozen_boundary_condition.hpp>
#include <zisa/boundary/no_boundary_condition.hpp>
#include <zisa/cli/input_parameters.hpp>

namespace zisa {

std::shared_ptr<BoundaryCondition>
make_boundary_condition(const InputParameters &params,
                        std::shared_ptr<Grid> grid,
                        std::shared_ptr<AllVariables> all_vars);

}
#endif // ZISA_BOUNDARY_CONDITION_FACTORY_HPP
