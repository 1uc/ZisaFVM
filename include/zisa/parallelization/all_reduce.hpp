#ifndef ZISA_ALL_REDUCE_HPP_CYQKU
#define ZISA_ALL_REDUCE_HPP_CYQKU

#include <zisa/config.hpp>
#include <zisa/parallelization/mpi.hpp>

namespace zisa {

enum class ReductionOperation { min };

class AllReduce {
public:
  virtual ~AllReduce() = default;

  double operator()(double local) const;

protected:
  virtual double do_reduce(double local) const = 0;
};

}
#endif // ZISA_ALL_REDUCE_HPP