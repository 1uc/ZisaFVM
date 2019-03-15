#include <fstream>

#include <zisa/testing/testing_framework.hpp>

#include <zisa/io/file_manipulation.hpp>
#include <zisa/io/file_name_generator.hpp>

TEST_CASE("FileNameGenerator; find_last_data_file", "[io]") {
  auto fng = zisa::FileNameGenerator("data/__fng", "_data-%04d", ".h5");

  zisa::create_directory(zisa::dirname(fng.next_name()));
  fng.advance_to(10);

  std::vector<std::string> filenames;
  filenames.reserve(7);
  for (int i = 0; i < 5; ++i) {
    filenames.emplace_back(fng.next_name());
  }
  filenames.emplace_back(fng.steady_state_filename);
  filenames.emplace_back(fng.grid_filename);

  for (const auto &filename : filenames) {
    std::ofstream os(filename);
    assert(os.good());

    os << ":)";
  }

  auto actual = zisa::find_last_data_file(fng);
  auto expected = std::string("data/__fng_data-0014.h5");

  REQUIRE(actual == expected);
}