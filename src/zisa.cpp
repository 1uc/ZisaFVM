/* Initialize libraries and launch the simulation.
 */

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <unistd.h>
#include <zisa/experiments/down_sample_reference.hpp>

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

void run_zisa(const std::string &mode, const zisa::InputParameters &params) {
  auto experiment = zisa::make_experiment(params);

  if (mode == "run") {
    experiment->run();
  } else if (mode == "post-process") {
    experiment->post_process();
  }
}

int main(int argc, char *argv[]) {
  signal(SIGSEGV, handler); // install our handler

#if ZISA_HAS_MPI == 1
  int requested = MPI_THREAD_MULTIPLE;
  int provided = -1;

  MPI_Init_thread(&argc, &argv, requested, &provided);
  LOG_ERR_IF(requested > provided,
             "MPI does not support the requested level of multi-threading.");
#endif

#if ZISA_HAS_OPENGL == 1
  glewExperimental = true;
  auto glfw_status = glfwInit();
  LOG_ERR_IF(!glfw_status, "Failed to initialize GLFW.\n");
#endif

  auto [mode, input_parameters] = zisa::parse_command_line(argc, argv);
  run_zisa(mode, input_parameters);

#if ZISA_HAS_MPI == 1
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
