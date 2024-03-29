// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef PARSE_COMMAND_LINE_H_WIQC2
#define PARSE_COMMAND_LINE_H_WIQC2

#include <zisa/cli/input_parameters.hpp>

namespace zisa {

std::pair<std::string, InputParameters> parse_command_line(int argc,
                                                           char *argv[]);

} // namespace zisa

#endif /* end of include guard */
