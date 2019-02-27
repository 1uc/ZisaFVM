#ifndef FILE_NAME_GENERATOR_H_VLEA4
#define FILE_NAME_GENERATOR_H_VLEA4

#include <string>

#include <zisa/config.hpp>

namespace zisa {

/// Sequentially numbered file names.
/** Provides consistent names for snapshots, the grid and the steady-state
 *  files. The snapshot names are a incrementally numbered sequence.
 *
 *  Example:
 *
 *      auto fng = FileNameGenerator("data/foo", "-%04d", ".h5");
 *
 *      std::cout << fng.grid_filename << "\n";
 *      std::cout << fng.steady_state_filename << "\n";
 *      for (int i = 0; i < 293; ++i) {
 *        std::cout << fng() << "\n";
 *      }
 *  will produce:
 *
 *      data/foo_grid.h5
 *      data/foo_steady-state.h5
 *      data/foo-0000.h5
 *        ...
 *      data/foo-0292.h5
 *
 *  @note The pattern should be a printf-style format for one integer.
 */
class FileNameGenerator {
public:
  FileNameGenerator(const std::string &stem,
                    const std::string &pattern,
                    const std::string &suffix);

  /// Generate the next numbered file name.
  std::string next_name(void);

  /// Generate numbers starting from `k`.
  void advance_to(int k);

  const std::string filename_stem;         ///< First part of all filenames.
  const std::string steady_state_filename; ///< Name of the steady-state.
  const std::string reference_filename;    ///< Name of the reference solution.
  const std::string grid_filename;         ///< Name of the grid.
  const std::string xdmf_grid_filename;    ///< Name of the grid.

private:
  std::string pattern;
  int count;
};

template <class Map>
std::shared_ptr<FileNameGenerator> make_file_name_generator(const Map &map) {
  return std::make_shared<FileNameGenerator>(
      map["stem"], map["pattern"], map["suffix"]);
}

} // namespace zisa
#endif /* end of include guard */
