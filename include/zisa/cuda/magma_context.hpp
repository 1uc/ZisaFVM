#ifndef MAGMA_CONTEXT_HPP
#define MAGMA_CONTEXT_HPP

#include <zisa/config.hpp>

#include <magma_v2.h>
#include <memory>

namespace zisa {
namespace cuda {
namespace magma {

/// RAII style MAGMA context.
/** MAGMA requires a `magma_init()` and `magma_finalize()`.
 */
class MAGMAContext {
public:
  MAGMAContext() { magma_init(); }

  ~MAGMAContext() { magma_finalize(); }

  MAGMAContext(const MAGMAContext &) = delete;
  MAGMAContext(MAGMAContext &&) = default;

  const MAGMAContext &operator=(const MAGMAContext &) = delete;
  MAGMAContext &operator=(MAGMAContext &&) = default;
};

/// RAII style MAGMA queue.
/**  A MAGMA queue is one of those pointer to opaque object constructs
 *  frequently used in C. We here wrap it such that it adheres to the RAII
 *  principle.
 *
 *  Therefore, this class is non-copyable. If there are multiple users of the
 *  same queue, consider a smart pointer.
 */
class MAGMAQueue {
public:
  MAGMAQueue() = default;
  MAGMAQueue(std::shared_ptr<MAGMAContext> magma_context,
             magma_int_t device_id);

  MAGMAQueue(const MAGMAQueue &) = delete;
  MAGMAQueue(MAGMAQueue &&) = default;

  ~MAGMAQueue();

  const MAGMAQueue &operator=(const MAGMAQueue &) = delete;
  MAGMAQueue &operator=(MAGMAQueue &&) = default;

  inline void sync() const { magma_queue_sync(queue()); }
  inline const magma_queue_t &queue() const { return queue_; }

private:
  std::shared_ptr<MAGMAContext> magma_context;
  magma_queue_t queue_ = nullptr;
};

std::shared_ptr<MAGMAQueue> make_default_queue();

}
}
}
#endif // MAGMA_CONTEXT_HPP
