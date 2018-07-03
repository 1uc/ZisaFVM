/* Initialize libraries and launch the simulation.
 */

#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <zisa/utils/logging.hpp>
#include <zisa/cli/parse_command_line.hpp>

// TODO refactor
#include <zisa/grid/grid.hpp>

void handler(int sig) {
  void *array[10];

  // get void*'s for all entries on the stack
  auto size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void run_zisa(const zisa::InputParameters &params) {

  std::string grid_file = params["grid"]["file"];
  auto mesh = zisa::load_gmsh(grid_file);

  return;
}

int main(int argc, char *argv[])
{
  signal(SIGSEGV, handler); // install our handler

  run_zisa(zisa::parse_command_line(argc, argv));
  return EXIT_SUCCESS;
}
