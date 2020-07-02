#include <zisa/model/gravity.hpp>

namespace zisa {

RadialAlignment::RadialAlignment(double epsilon) : epsilon(epsilon) {}

void save(HDF5Writer &writer, const RadialAlignment &alignment) {
  writer.write_scalar(alignment.epsilon, "epsilon");
}

RadialAlignment RadialAlignment::load(HDF5Reader &reader) {
  auto eps = reader.read_scalar<double>("epsilon");
  return RadialAlignment(eps);
}

void save(HDF5Writer &writer, const AxialAlignment &alignment) {
  writer.write_scalar(alignment.axis[0], "x");
  writer.write_scalar(alignment.axis[1], "y");
  writer.write_scalar(alignment.axis[2], "z");
}

AxialAlignment AxialAlignment::load(HDF5Reader &reader) {
  auto x = reader.read_scalar<double>("x");
  auto y = reader.read_scalar<double>("y");
  auto z = reader.read_scalar<double>("z");

  return AxialAlignment(XYZ{x, y, z});
}

void save(HDF5Writer &writer, const ConstantGravity &gravity) {
  writer.write_scalar(gravity.gravity, "g");
}

ConstantGravity ConstantGravity::load(HDF5Reader &reader) {
  auto g_ = reader.read_scalar<double>("g");
  return ConstantGravity(g_);
}

PointMassGravity::PointMassGravity(double GM, double X) : GM(GM), X(X) {}

PointMassGravity::PointMassGravity(double G, double M, double X)
    : GM(G * M), X(X) {}

void save(HDF5Writer &writer, const PointMassGravity &gravity) {
  writer.write_scalar(gravity.GM, "GM");
  writer.write_scalar(gravity.X, "X");
}

PointMassGravity PointMassGravity::load(HDF5Reader &reader) {
  auto GM_ = reader.read_scalar<double>("GM");
  auto X_ = reader.read_scalar<double>("X");

  return PointMassGravity(GM_, X_);
}

void save(HDF5Writer &writer, const PolytropeGravity &gravity) {
  writer.write_scalar(gravity.rhoC, "rhoC");
  writer.write_scalar(gravity.K, "K");
  writer.write_scalar(gravity.G, "G");
  writer.write_scalar(gravity.eps, "eps");
}

PolytropeGravity::PolytropeGravity(double rhoC, double K, double G)
    : rhoC(rhoC), K(K), G(G) {}

PolytropeGravity::PolytropeGravity(double rhoC, double K, double G, double eps)
    : rhoC(rhoC), K(K), G(G), eps(eps) {}

PolytropeGravity PolytropeGravity::load(HDF5Reader &reader) {
  auto rhoC_ = reader.read_scalar<double>("rhoC");
  auto K_ = reader.read_scalar<double>("K");
  auto G_ = reader.read_scalar<double>("G");
  auto eps_ = reader.read_scalar<double>("eps");

  return PolytropeGravity(rhoC_, K_, G_, eps_);
}

void save(HDF5Writer &writer, const PolytropeGravityWithJump &gravity) {
  writer.write_scalar(gravity.r_crit, "r_crit");

  writer.open_group("inner");
  save(writer, gravity.inner);
  writer.switch_group("outer");
  save(writer, gravity.outer);
  writer.close_group();
}

SphericalGravity::SphericalGravity(array<double, 1> radii) {
  auto phi = array<double, 1>(radii.shape());
  interpolate
      = NonUniformLinearInterpolation<double>(std::move(radii), std::move(phi));
}

SphericalGravity::SphericalGravity(array<double, 1> radii, array<double, 1> phi)
    : interpolate(std::move(radii), std::move(phi)) {}

void save(HDF5Writer &writer, const SphericalGravity &gravity) {
  LOG_WARN("Not saving the gravitational potential.");
  //  save(writer, gravity.interpolate.points, "radii");
  //  save(writer, gravity.interpolate.values, "phi");
}

SphericalGravity SphericalGravity::load(HDF5Reader &reader) {
  auto radii_ = array<double, 1>::load(reader, "radii");
  auto phi_ = array<double, 1>::load(reader, "phi");

  return SphericalGravity(std::move(radii_), std::move(phi_));
}

void save(HDF5Writer &, const NoGravity &) { return; }

} // zisa
