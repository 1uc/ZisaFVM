// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#ifndef FILE_NAME_GENERATOR_H_VLEA4
#define FILE_NAME_GENERATOR_H_VLEA4

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include <zisa/config.hpp>
#include <zisa/io/file_manipulation.hpp>

namespace zisa {

/// Minimal interface of a file name generator.
class FNG {
public:
  virtual std::string next_name() = 0;
  virtual std::string steady_state() = 0;
};

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
 *      data/grid.h5
 *      data/steady-state.h5
 *      data/foo-0000.h5
 *        ...
 *      data/foo-0292.h5
 *
 *  @note The pattern should be a printf-style format for one integer.
 */
class FileNameGenerator : public FNG {
public:
  FileNameGenerator(const std::string &stem,
                    const std::string &pattern,
                    const std::string &suffix);

  FileNameGenerator(const std::string &dir,
                    const std::string &stem,
                    const std::string &pattern,
                    const std::string &suffix);

  /// Generate the filename of a specified generation.
  std::string filename(int generation);

  /// Generate the next numbered file name.
  virtual std::string next_name() override;

  /// Generate the file name for the steady state.
  virtual std::string steady_state() override;

  /// Generate numbers starting from `k`.
  void advance_to(int k);
  void advance_to(const std::string &filename);

  /// Returns the count/generation of a datafile.
  /** Returns -1 if it's not a datafile.
   */
  int generation(const std::filesystem::path &path);

  const std::string filename_stem;         ///< First part of all filenames.
  const std::string steady_state_filename; ///< Path of the steady-state.
  const std::string grid_filename;         ///< Path of the grid.
  const std::string dirname;               ///< Directory containing all files.

private:
  std::string pattern_;
  int count_;
};

class SingleFileNameGenerator : public FNG {
public:
  SingleFileNameGenerator(std::string filename)
      : filename(std::move(filename)) {}

  virtual std::string next_name() override {
    LOG_ERR_IF(used, "This filename generator has used up all its names.");
    return filename;
  }

  virtual std::string steady_state() override {
    LOG_ERR("No steady state filename available.");
  }

private:
  std::string filename;
  bool used = false;
};

class FixedFileNameGenerator : public FNG {
public:
  FixedFileNameGenerator(std::vector<std::string> filenames)
      : filenames(std::move(filenames)) {}

  virtual std::string next_name() override {
    LOG_ERR_IF(n_used >= filenames.size(),
               "This filename generator has used up all its names.");

    return filenames[n_used];
  }

  virtual std::string steady_state() override {
    LOG_ERR("No steady state filename available.");
  }

private:
  std::vector<std::string> filenames;
  std::size_t n_used = 0;
};

std::shared_ptr<FileNameGenerator>
make_file_name_generator(const std::string &dir,
                         const std::string &stem,
                         const std::string &pattern,
                         const std::string &suffix);

std::string find_last_data_file(FileNameGenerator &fng);
std::string find_first_data_file(FileNameGenerator &fng);

} // namespace zisa
#endif /* end of include guard */
