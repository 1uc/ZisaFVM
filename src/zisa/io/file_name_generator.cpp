#include <zisa/io/file_name_generator.hpp>

namespace zisa {
FileNameGenerator::FileNameGenerator(const std::string &stem,
                                     const std::string &pattern,
                                     const std::string &suffix)
    : FileNameGenerator("", stem, pattern, suffix) {}

FileNameGenerator::FileNameGenerator(const std::string &dir,
                                     const std::string &stem,
                                     const std::string &pattern,
                                     const std::string &suffix)
    : filename_stem(stem),
      steady_state_filename(dir + "steady_state" + suffix),
      reference_filename(dir + stem + "_reference" + suffix),
      grid_filename(dir + "grid" + suffix),
      pattern_(dir + stem + pattern + suffix),
      count_(0) {}

std::string FileNameGenerator::filename(int generation) {
  return string_format(pattern_, generation);
}


std::string FileNameGenerator::next_name() {
  return filename(count_++);
}

void FileNameGenerator::advance_to(int k) { count_ = k; }
void FileNameGenerator::advance_to(const std::string &filename) {
  int k = -1;
  sscanf(filename.c_str(), pattern_.c_str(), &k);

  advance_to(k + 1);
}

int FileNameGenerator::generation(const std::filesystem::path &rel_path) {
  int gen = -1;

  auto path = std::string(std::filesystem::absolute(rel_path));
  auto pattern = std::string(std::filesystem::absolute(pattern_));
  auto status = sscanf(path.c_str(), pattern.c_str(), &gen);

  return (status == 1 ? gen : -1);
}

std::string find_last_data_file(FileNameGenerator &fng) {
  namespace fs = std::filesystem;

  std::string path = fs::relative(zisa::dirname(fs::absolute(fng.filename_stem)));
  LOG_ERR_IF(!fs::is_directory(path),
             string_format("Not a directory. [%s]", path.c_str()));

  auto files = std::vector<fs::path>{};
  for (const auto &entry : fs::directory_iterator(path)) {
    files.push_back(entry.path());
  }

  auto m = std::max_element(
      begin(files), end(files), [&fng](const fs::path &p1, const fs::path &p2) {
        return fng.generation(p1) < fng.generation(p2);
      });

  LOG_ERR_IF(fng.generation(*m) < 0, "No data-files.");
  return std::string(fs::relative(*m));
};

std::string find_first_data_file(FileNameGenerator &fng) {
  return fng.filename(0);
};


std::shared_ptr<FileNameGenerator>
make_file_name_generator(const std::string &dir,
                         const std::string &stem,
                         const std::string &pattern,
                         const std::string &suffix) {
  return std::make_shared<FileNameGenerator>(dir, stem, pattern, suffix);
}

} // namespace zisa
