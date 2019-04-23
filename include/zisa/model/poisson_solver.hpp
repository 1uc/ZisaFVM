#ifndef ZISA_POISSON_SOLVER_HPP
#define ZISA_POISSON_SOLVER_HPP

#include <zisa/model/all_variables_fwd.hpp>
#include <zisa/model/gravity.hpp>

namespace zisa {

class PoissonSolver {
public:
  virtual ~PoissonSolver() = default;
  virtual void update(RadialGravity &gravity,
                      const AllVariables &all_variables) const = 0;
};

}

#endif // ZISA_POISSON_SOLVER_HPP
