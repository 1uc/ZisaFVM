#include <zisa/model/gravity.hpp>

namespace zisa {

SphericalGravity::SphericalGravity(array<double, 1> radii) {
  auto phi = array<double, 1>(radii.shape());
  interpolate = std::make_shared<NonUniformLinearInterpolation<double>>(
      std::move(radii), std::move(phi));
}

SphericalGravity::SphericalGravity(array<double, 1> radii, array<double, 1> phi)
    : interpolate(std::make_shared<NonUniformLinearInterpolation<double>>(
          std::move(radii), std::move(phi))) {}

void save(HDF5Writer &writer, const RadialAlignment &alignment) {
  writer.write_scalar(alignment.epsilon, "epsilon");
}

void save(HDF5Writer &writer, const AxialAlignment &alignment) {
  writer.write_scalar(alignment.axis[0], "x");
  writer.write_scalar(alignment.axis[1], "y");
  writer.write_scalar(alignment.axis[2], "z");
}

void save(HDF5Writer &writer, const ConstantGravity &gravity) {
  writer.write_scalar(gravity.gravity, "g");
}

void save(HDF5Writer &writer, const PointMassGravity &gravity) {
  writer.write_scalar(gravity.GM, "GM");
  writer.write_scalar(gravity.X, "X");
}

void save(HDF5Writer &writer, const PolytropeGravity &gravity) {
  writer.write_scalar(gravity.rhoC, "rhoC");
  writer.write_scalar(gravity.K, "K");
  writer.write_scalar(gravity.G, "G");
  writer.write_scalar(gravity.eps, "eps");
}

void save(HDF5Writer &writer, const PolytropeGravityWithJump &gravity) {
  writer.write_scalar(gravity.r_crit, "r_crit");

  writer.open_group("inner");
  save(writer, gravity.inner);
  writer.switch_group("outer");
  save(writer, gravity.outer);
  writer.close_group();
}

void save(HDF5Writer &writer, const SphericalGravity &gravity) {
  save(writer, gravity.interpolate->points, "radii");
  save(writer, gravity.interpolate->values, "phi");
}

} // zisa
