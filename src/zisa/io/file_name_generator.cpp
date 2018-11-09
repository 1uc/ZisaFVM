#include <zisa/io/file_name_generator.hpp>

namespace zisa {

FileNameGenerator::FileNameGenerator(const std::string &stem,
                                     const std::string &pattern,
                                     const std::string &suffix)
    : filename_stem(stem),
      steady_state_filename(stem + "_steady-state" + suffix),
      reference_filename(stem + "_reference" + suffix),
      grid_filename(stem + "_grid" + suffix),
      xdmf_grid_filename(stem + "_xdmf_grid" + suffix),
      pattern(stem + pattern + suffix),
      count(0) {}

std::string FileNameGenerator::next_name(void) {
  std::string file_name = string_format(pattern, count);

  ++count;
  return file_name;
}

void FileNameGenerator::advance_to(int k) { count = k; }

} // namespace zisa
