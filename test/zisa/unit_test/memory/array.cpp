#include <zisa/io/hdf5_serial_writer.hpp>
#include <zisa/memory/array.hpp>
#include <zisa/memory/array_view.hpp>
#include <zisa/memory/column_major.hpp>
#include <zisa/testing/testing_framework.hpp>

using namespace zisa;

bool implicit_conversion(const array_const_view<double, 3> &view) {

  return view.raw() != nullptr;
}

TEST_CASE("array; basics", "[array]") {
  auto a = array<double, 3>({3ul, 3ul, 2ul}, device_type::cpu);
  auto b = array<double, 3>({3ul, 3ul, 2ul}, device_type::cpu);

  REQUIRE(a.raw() != static_cast<const decltype(b) &>(b).raw());
  REQUIRE(a.shape() == b.shape());

  auto b_view = array_view<double, 3>(b);
  auto bb_view = array_view<double, 3>(shape(b), raw_ptr(b));
  auto bc_view = array_const_view<double, 3>(shape(b), raw_ptr(b));

  REQUIRE(b.raw() == b_view.raw());

  b(1, 0, 1) = -42.0;
  REQUIRE(b_view(1, 0, 1) == -42.0);

  b_view(0, 0, 0) = 42.0;
  REQUIRE(b(0, 0, 0) == 42.0);

  REQUIRE(implicit_conversion(a));
}

TEST_CASE("array; write to file", "[array]") {

  auto filename = "__unit_tests--array-to-hdf5.h5";
  auto label = "a";

  auto shape = zisa::shape_t<3>{3ul, 4ul, 2ul};

  auto a = array<double, 3>(shape);
  for (zisa::int_t i = 0; i < a.size(); ++i) {
    a[i] = i;
  }

  {
    auto writer = HDF5SerialWriter(filename);
    zisa::save(writer, a, label);
  }

  {
    auto reader = HDF5SerialReader(filename);
    auto dims = reader.dims(label);

    for (zisa::int_t k = 0; k < 3; ++k) {
      REQUIRE(shape[0] == static_cast<zisa::int_t>(dims[0]));
    }

    auto b = zisa::array<double, 3>::load(reader, label);
    REQUIRE(b == a);
  }
}

TEST_CASE("array; builds for general Indexing.", "[array]") {

  // The point is to check that `array<double, 3, Indexing>`
  // compiles fine. Despite the fact that `save` and `load` only
  // work if `Indexing == row_major`.

  auto shape = zisa::shape_t<3>{3ul, 4ul, 2ul};

  auto a = array<double, 3, column_major>(shape);
  for (zisa::int_t i = 0; i < a.size(); ++i) {
    a[i] = i;
  }
}
