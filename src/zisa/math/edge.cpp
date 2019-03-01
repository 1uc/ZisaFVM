#include <zisa/math/edge.hpp>

namespace zisa {

std::ostream &operator<<(std::ostream &os, const Edge &edge) {
  os << "Edge(" << edge.start_point() << ", " << edge.end_point() << ")";
  return os;
}

}
