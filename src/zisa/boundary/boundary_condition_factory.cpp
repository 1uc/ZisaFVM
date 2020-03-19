#include <zisa/boundary/boundary_condition_factory.hpp>

namespace zisa {
std::shared_ptr<BoundaryCondition>
make_boundary_condition(const InputParameters &params,
                        std::shared_ptr<Grid> grid,
                        std::shared_ptr<AllVariables> all_vars) {

  if (has_key(params, "boundary-condition")) {
    const auto &bc = params["boundary-condition"].value("mode", "");
    if (bc == "frozen"){
      return std::make_shared<FrozenBC>(*grid, *all_vars);
    }
  }

  return std::make_shared<NoBoundaryCondition>();
}
}