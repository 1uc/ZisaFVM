#ifndef ZISA_HALO_EXCHANGE_BC_HPP_ICSQO
#define ZISA_HALO_EXCHANGE_BC_HPP_ICSQO

#include <zisa/boundary/boundary_condition.hpp>
#include <zisa/model/all_variables.hpp>
#include <zisa/parallelization/halo_exchange.hpp>
#include <zisa/utils/indent_block.hpp>

namespace zisa {

class HaloExchangeBC : public BoundaryCondition {
public:
  HaloExchangeBC(std::shared_ptr<BoundaryCondition> local_bc,
                 std::shared_ptr<HaloExchange> halo_exchange);

  void apply(AllVariables &u, double t) override;

  std::string str() const override;

private:
  std::shared_ptr<BoundaryCondition> local_bc;
  std::shared_ptr<HaloExchange> halo_exchange;
};

}

#endif // ZISA_HALO_EXCHANGE_BC_HPP
