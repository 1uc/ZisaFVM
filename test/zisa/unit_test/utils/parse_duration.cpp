
#include <chrono>
#include <vector>

#include <zisa/testing/testing_framework.hpp>
#include <zisa/utils/parse_duration.hpp>

TEST_CASE("parse_duration", "[utils]") {

  using namespace std::chrono_literals;
  auto examples
      = std::vector<std::pair<std::string, std::chrono::milliseconds>>{
          {"43s", 43s}, {"400ms", 400ms}, {"10ms", 10ms}};

  for (auto &[s, d] : examples) {
    INFO(s);
    REQUIRE(zisa::parse_duration_ms(s) == d);
  }
}
