/* Parse the command line arguments.
 */

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <zisa/cli/input_parameters.hpp>
#include <zisa/utils/logging.hpp>

namespace zisa {

InputParameters parse_command_line(int argc, char *argv[]) {
  po::variables_map options;

  // clang-format off
  // generic options
  po::options_description generic("Generic options");
  generic.add_options()("help,h", "produce this message")
    // ("version", "output version information")
    ("config,c", po::value<std::string>(), "config file to use")
    ;
  // clang-format on

  // command line ONLY options
  po::options_description cmdline("Command line options");
  cmdline.add(generic);

  // first parse cmdline and check what config file to use
  po::store(po::parse_command_line(argc, argv, cmdline), options);
  po::notify(options);

  LOG_ERR_IF(options.count("config") != 1, "No config file provided.");

  std::string config_file(options["config"].as<std::string>());
  return InputParameters(config_file);
}

} // namespace zisa
