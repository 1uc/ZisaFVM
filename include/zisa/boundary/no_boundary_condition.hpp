#ifndef NO_BOUNDARY_CONDITION_H_GKGJL
#define NO_BOUNDARY_CONDITION_H_GKGJL

#include <zisa/config.hpp>
#include <zisa/boundary/boundary_condition.hpp>

namespace zisa {

class NoBoundaryCondition : public BoundaryCondition {
public:
  virtual void apply(AllVariables &u, double t) override;
  virtual std::string str() const override;
};

} // namespace zisa

#endif /* end of include guard */
