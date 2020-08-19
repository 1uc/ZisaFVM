#include <zisa/testing/testing_framework.hpp>

#include <filesystem>

TEST_CASE("std::filesystem; join", "[fs]") {
  auto a = std::string("");
  auto b = std::string("foo");

  REQUIRE(std::filesystem(a)/b == std::filesystem(b));
}
