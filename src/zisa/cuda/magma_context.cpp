#include <zisa/cuda/magma_context.hpp>

namespace zisa {
namespace cuda {
namespace magma {

MAGMAQueue::MAGMAQueue(std::shared_ptr<MAGMAContext> magma_context_,
                       magma_int_t device_id)
    : magma_context(std::move(magma_context_)) {
  magma_queue_create(device_id, &queue_);
}

MAGMAQueue::~MAGMAQueue() {
  if (queue_ != nullptr) {
    magma_queue_destroy(queue_);
  }
}

std::shared_ptr<MAGMAQueue> make_default_queue() {
  return std::make_shared<MAGMAQueue>(std::make_shared<MAGMAContext>(), 0);
}
}
}
}
