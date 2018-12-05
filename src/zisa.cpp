/* Initialize libraries and launch the simulation.
 */

#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#if ZISA_HAS_OPENGL == 1
// -- must come before all other GL stuff.
#include <GL/glew.h>
// ---------------------------------------

#include <GLFW/glfw3.h>
#endif

#include <zisa/cli/parse_command_line.hpp>
#include <zisa/utils/logging.hpp>

#include <zisa/experiments/numerical_experiment_factory.hpp>

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

  auto experiment = zisa::make_experiment(params);
  experiment->run();
}

int main(int argc, char *argv[]) {
  signal(SIGSEGV, handler); // install our handler

#if ZISA_HAS_OPENGL == 1
  glewExperimental = true;
  if (!glfwInit()) {
    std::cerr << "Failed to initialize GLFW.\n";
    return -1;
  }
#endif

  run_zisa(zisa::parse_command_line(argc, argv));
  return EXIT_SUCCESS;
}
