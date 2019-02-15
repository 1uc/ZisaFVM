#include <zisa/testing/testing_framework.hpp>

#include <zisa/memory/array.hpp>
#include <zisa/memory/block_allocator.hpp>

TEST_CASE("block_allocator; basic API", "[memory]") {

  auto shape = zisa::shape_t<2>{3ul, 4ul};

  std::size_t max_elements = 10;
  auto allocator
      = std::make_shared<zisa::block_allocator<zisa::array<double, 2>>>(
          max_elements);

  SECTION("Check reuse.") {
    void *address_call_1;
    {
      auto ptr = allocator->allocate(shape);
      address_call_1 = &(*ptr);
    }

    void *address_call_2;
    {
      auto ptr = allocator->allocate(shape);
      address_call_2 = &(*ptr);
    }

    REQUIRE(address_call_2 == address_call_1);
  }

  SECTION("Check distinct.") {
    auto ptr1 = allocator->allocate(shape);
    auto ptr2 = allocator->allocate(shape);

    REQUIRE(&(*ptr1) != &(*ptr2));
  }

  SECTION("Throw when out of memory.") {
    std::vector<zisa::locked_ptr<zisa::array<double, 2>>> ptrs;

    for (std::size_t i = 0; i < max_elements; ++i) {
      ptrs.push_back(allocator->allocate(shape));
    }

    REQUIRE_THROWS(ptrs.push_back(allocator->allocate(shape)));
  }

  SECTION("Call frequently.") {
    for (std::size_t i = 0; i < 10 * max_elements; ++i) {
      auto ptr = allocator->allocate(shape);
      REQUIRE(ptr->raw() != nullptr);
    }
  }
}
