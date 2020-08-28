#include <zisa/testing/testing_framework.hpp>

#include <filesystem>

TEST_CASE("std::filesystem; join", "[fs]") {
  auto a = std::string("");
  auto b = std::string("foo");

  REQUIRE(std::filesystem::path(a) / b == std::filesystem::path(b));
}
