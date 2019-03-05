/* Parse the command line arguments.
 */

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <zisa/cli/input_parameters.hpp>
#include <zisa/utils/logging.hpp>

namespace zisa {

std::pair<std::string, InputParameters> parse_command_line(int argc,
                                                           char *argv[]) {

  // Description & thanks to: https://stackoverflow.com/a/23098581

  po::variables_map options;

  // clang-format off
  // generic options
  po::options_description generic("Generic options");
  generic.add_options()
      ("help,h", "produce this message")
      ("subcmd", "subcommand to run")
      ("subargs", "arguments for the subcommand")
    // ("version", "output version information")
    ;
  po::positional_options_description pos;
  pos.add("subcmd", 1).add("subargs", -1);
  // clang-format on

  // first parse cmdline and check what config file to use
  po::parsed_options parsed = po::command_line_parser(argc, argv)
                                  .options(generic)
                                  .positional(pos)
                                  .allow_unregistered()
                                  .run();

  po::store(parsed, options);

  if (options.count("help") != 0) {
    std::cout << generic << "\n";
    std::exit(EXIT_SUCCESS);
  }

  po::options_description config("Options");
  config.add_options()(
      "config,c", po::value<std::string>(), "config file to use");

  LOG_ERR_IF(options.count("subcmd") != 1, "No subcommand provided.");
  std::string subcmd = options["subcmd"].as<std::string>();

  std::vector<std::string> sub_opts
      = po::collect_unrecognized(parsed.options, po::include_positional);
  sub_opts.erase(sub_opts.begin());

  po::store(po::command_line_parser(sub_opts).options(config).run(), options);
  po::notify(options);

  LOG_ERR_IF(options.count("config") != 1, "No config file provided.");
  std::string config_file(options["config"].as<std::string>());

  return {subcmd, InputParameters(config_file)};
}

} // namespace zisa
