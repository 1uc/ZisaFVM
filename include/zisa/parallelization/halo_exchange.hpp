// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef ZISA_HALO_EXCHANGE_HPP_IDOIW
#define ZISA_HALO_EXCHANGE_HPP_IDOIW

#include <zisa/config.hpp>
#include <zisa/model/all_variables.hpp>

namespace zisa {
class HaloExchange {
public:
  virtual ~HaloExchange() = default;

  /// Post the halo exchange.
  /** The exchange is not guaranteed to be complete until `wait` is called. */
  virtual void operator()(AllVariables &all_vars) = 0;

  /// Wait for the halo exchange to be completed.
  virtual void wait() = 0;
};

class NoHaloExchange : public HaloExchange {
public:
  virtual ~NoHaloExchange() = default;

  virtual void operator()(AllVariables &all_vars) override;
  virtual void wait() override;
};

}
#endif // ZISA_HALO_EXCHANGE_HPP
