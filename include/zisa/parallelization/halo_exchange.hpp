#ifndef ZISA_HALO_EXCHANGE_HPP_IDOIW
#define ZISA_HALO_EXCHANGE_HPP_IDOIW

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {
class HaloExchange {
public:
  virtual ~HaloExchange() = default;
  virtual void operator()(AllVariables &all_vars) = 0;
};

}
#endif // ZISA_HALO_EXCHANGE_HPP
