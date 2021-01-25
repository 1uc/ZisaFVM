#include <zisa/parallelization/halo_exchange.hpp>

namespace zisa {
void NoHaloExchange::operator()(AllVariables &) { return; }
void NoHaloExchange::wait() { return; }
}
