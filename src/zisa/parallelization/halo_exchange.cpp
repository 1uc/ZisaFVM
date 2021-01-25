#include <zisa/parallelization/halo_exchange.hpp>

namespace zisa {
void NoHaloExchange::operator()(AllVariables &all_vars) { return; }
void NoHaloExchange::wait() { return; }
}
