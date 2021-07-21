// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef BACKTRACE_H_DMH02QWB
#define BACKTRACE_H_DMH02QWB

#include <string>

namespace zisa {
/// Print a backtrace with line numbers.
/** Thanks:
 *     http://stackoverflow.com/a/4611112/5103043
 */
std::string backtrace(void);

}

#endif /* end of include guard: BACKTRACE_H_DMH02QWB */
