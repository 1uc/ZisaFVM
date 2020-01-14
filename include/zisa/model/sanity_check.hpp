/* Check that the array contains sane values.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-09-06
 */
#ifndef SANITY_CHECK_H_L4SIPN6H
#define SANITY_CHECK_H_L4SIPN6H

#include "zisa/config.hpp"
#include "zisa/model/all_variables.hpp"

namespace zisa {

class SanityCheck {
public:
  virtual ~SanityCheck() = default;
  virtual bool operator()(const AllVariables &all_variables) const = 0;
};

class NoSanityCheck : public SanityCheck {
public:
  bool operator()(const AllVariables &) const override {
    return true;
  }
};

} // namespace zisa

#endif /* end of include guard: SANITY_CHECK_H_LPN4SI6H */
