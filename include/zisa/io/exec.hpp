#ifndef EXEC_H_GIT5YS0K
#define EXEC_H_GIT5YS0K
#include <string>

namespace zisa {

/// Execute a command in Linux.
/** Thanks:
 *    http://stackoverflow.com/a/478960/5103043
 */
std::string exec(const std::string &cmd);

}

#endif /* end of include guard: EXEC_H_GIT5YS0K */
