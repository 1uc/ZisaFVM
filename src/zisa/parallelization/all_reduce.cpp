#include <zisa/parallelization/all_reduce.hpp>

namespace zisa {

double AllReduce::operator()(double local) const { return do_reduce(local); }

}