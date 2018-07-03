/* Simple Linux backtrace with line numbers.
 *
 * Authors: Luc Grosheintz <forbugrep@zoho.com>
 *    Date: 2016-01-14
 */
#ifndef BACKTRACE_H_DMH02QWB
#define BACKTRACE_H_DMH02QWB

#include <string>

namespace zisa {
/// Print a backtrace with line numbers.
/** Thanks:
 *     http://stackoverflow.com/a/4611112/5103043
 */
std::string backtrace(void);

} // namespace tyr

#endif /* end of include guard: BACKTRACE_H_DMH02QWB */
