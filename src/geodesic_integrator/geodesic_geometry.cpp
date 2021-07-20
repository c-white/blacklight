// Blacklight geodesic integrator - coordinates and metric components

// C++ headers
#include <cmath>  // hypot, sqrt

// Blacklight headers
#include "geodesic_integrator.hpp"
#include "../utils/array.hpp"       // Array

//--------------------------------------------------------------------------------------------------

// Function for calculating radial coordinate given location in coordinates used for geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   returned value: radial coordinate
// Notes:
//   Assumes Cartesian Kerr-Schild coordinates.
double GeodesicIntegrator::RadialGeodesicCoordinate(double x, double y, double z)
{
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  return r;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating covariant metric components for integrating geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   gcov: components set
// Notes:
//   Assumes gcov is allocated to be 4*4.
//   Assumes Cartesian Kerr-Schild coordinates (assumes Minkowski coordinates if ray_flat == true).
void GeodesicIntegrator::CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov)
{
  // Handle flat case
  if (ray_flat)
  {
    gcov.Zero();
    gcov(0,0) = -1.0;
    gcov(1,1) = 1.0;
    gcov(2,2) = 1.0;
    gcov(3,3) = 1.0;
    return;
  }

  // Calculate useful quantities
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r2 * r / (r2 * r2 + a2 * z * z);

  // Calculate null vector
  double l_0 = 1.0;
  double l_1 = (r * x + bh_a * y) / (r2 + a2);
  double l_2 = (r * y - bh_a * x) / (r2 + a2);
  double l_3 = z / r;

  // Calculate metric components
  gcov(0,0) = f * l_0 * l_0 - 1.0;
  gcov(0,1) = f * l_0 * l_1;
  gcov(0,2) = f * l_0 * l_2;
  gcov(0,3) = f * l_0 * l_3;
  gcov(1,0) = f * l_1 * l_0;
  gcov(1,1) = f * l_1 * l_1 + 1.0;
  gcov(1,2) = f * l_1 * l_2;
  gcov(1,3) = f * l_1 * l_3;
  gcov(2,0) = f * l_2 * l_0;
  gcov(2,1) = f * l_2 * l_1;
  gcov(2,2) = f * l_2 * l_2 + 1.0;
  gcov(2,3) = f * l_2 * l_3;
  gcov(3,0) = f * l_3 * l_0;
  gcov(3,1) = f * l_3 * l_1;
  gcov(3,2) = f * l_3 * l_2;
  gcov(3,3) = f * l_3 * l_3 + 1.0;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating contravariant metric components for integrating geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   gcon: components set
// Notes:
//   Assumes gcon is allocated to be 4*4.
//   Assumes Cartesian Kerr-Schild coordinates (assumes Minkowski coordinates if ray_flat == true).
void GeodesicIntegrator::ContravariantGeodesicMetric(double x, double y, double z,
    Array<double> &gcon)
{
  // Handle flat case
  if (ray_flat)
  {
    gcon.Zero();
    gcon(0,0) = -1.0;
    gcon(1,1) = 1.0;
    gcon(2,2) = 1.0;
    gcon(3,3) = 1.0;
    return;
  }

  // Calculate useful quantities
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r2 * r / (r2 * r2 + a2 * z * z);

  // Calculate null vector
  double l0 = -1.0;
  double l1 = (r * x + bh_a * y) / (r2 + a2);
  double l2 = (r * y - bh_a * x) / (r2 + a2);
  double l3 = z / r;

  // Calculate metric components
  gcon(0,0) = -f * l0 * l0 - 1.0;
  gcon(0,1) = -f * l0 * l1;
  gcon(0,2) = -f * l0 * l2;
  gcon(0,3) = -f * l0 * l3;
  gcon(1,0) = -f * l1 * l0;
  gcon(1,1) = -f * l1 * l1 + 1.0;
  gcon(1,2) = -f * l1 * l2;
  gcon(1,3) = -f * l1 * l3;
  gcon(2,0) = -f * l2 * l0;
  gcon(2,1) = -f * l2 * l1;
  gcon(2,2) = -f * l2 * l2 + 1.0;
  gcon(2,3) = -f * l2 * l3;
  gcon(3,0) = -f * l3 * l0;
  gcon(3,1) = -f * l3 * l1;
  gcon(3,2) = -f * l3 * l2;
  gcon(3,3) = -f * l3 * l3 + 1.0;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating derivatives of contravariant metric components for integrating geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   dgcon: components set
// Notes:
//   Assumes dgcon is allocated to be 3*4*4.
//   Assumes Cartesian Kerr-Schild coordinates (assumes Minkowski coordinates if ray_flat == true).
void GeodesicIntegrator::ContravariantGeodesicMetricDerivative(double x, double y, double z,
    Array<double> &dgcon)
{
  // Handle flat case
  if (ray_flat)
  {
    dgcon.Zero();
    return;
  }

  // Calculate useful quantities
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r2 * r / (r2 * r2 + a2 * z * z);

  // Calculate null vector
  double l0 = -1.0;
  double l1 = (r * x + bh_a * y) / (r2 + a2);
  double l2 = (r * y - bh_a * x) / (r2 + a2);
  double l3 = z / r;

  // Calculate scalar derivatives
  double dr_dx = r * x / (2.0 * r2 - rr2 + a2);
  double dr_dy = r * y / (2.0 * r2 - rr2 + a2);
  double dr_dz = (r * z + a2 * z / r) / (2.0 * r2 - rr2 + a2);
  double df_dx = -(r2 * r2 - 3.0 * a2 * z * z) * dr_dx / (r * (r2 * r2 + a2 * z * z)) * f;
  double df_dy = -(r2 * r2 - 3.0 * a2 * z * z) * dr_dy / (r * (r2 * r2 + a2 * z * z)) * f;
  double df_dz =
      -((r2 * r2 - 3.0 * a2 * z * z) * dr_dz + 2.0 * a2 * r * z) / (r * (r2 * r2 + a2 * z * z)) * f;

  // Calculate vector derivatives
  double dl0_dx = 0.0;
  double dl0_dy = 0.0;
  double dl0_dz = 0.0;
  double dl1_dx = ((x - 2.0 * r * l1) * dr_dx + r) / (r2 + a2);
  double dl1_dy = ((x - 2.0 * r * l1) * dr_dy + bh_a) / (r2 + a2);
  double dl1_dz = (x - 2.0 * r * l1) * dr_dz / (r2 + a2);
  double dl2_dx = ((y - 2.0 * r * l2) * dr_dx - bh_a) / (r2 + a2);
  double dl2_dy = ((y - 2.0 * r * l2) * dr_dy + r) / (r2 + a2);
  double dl2_dz = (y - 2.0 * r * l2) * dr_dz / (r2 + a2);
  double dl3_dx = -z / r2 * dr_dx;
  double dl3_dy = -z / r2 * dr_dy;
  double dl3_dz = -z / r2 * dr_dz + 1.0 / r;

  // Calculate metric component x-derivatives
  dgcon(0,0,0) = -(df_dx * l0 * l0 + f * dl0_dx * l0 + f * l0 * dl0_dx);
  dgcon(0,0,1) = -(df_dx * l0 * l1 + f * dl0_dx * l1 + f * l0 * dl1_dx);
  dgcon(0,0,2) = -(df_dx * l0 * l2 + f * dl0_dx * l2 + f * l0 * dl2_dx);
  dgcon(0,0,3) = -(df_dx * l0 * l3 + f * dl0_dx * l3 + f * l0 * dl3_dx);
  dgcon(0,1,0) = -(df_dx * l1 * l0 + f * dl1_dx * l0 + f * l1 * dl0_dx);
  dgcon(0,1,1) = -(df_dx * l1 * l1 + f * dl1_dx * l1 + f * l1 * dl1_dx);
  dgcon(0,1,2) = -(df_dx * l1 * l2 + f * dl1_dx * l2 + f * l1 * dl2_dx);
  dgcon(0,1,3) = -(df_dx * l1 * l3 + f * dl1_dx * l3 + f * l1 * dl3_dx);
  dgcon(0,2,0) = -(df_dx * l2 * l0 + f * dl2_dx * l0 + f * l2 * dl0_dx);
  dgcon(0,2,1) = -(df_dx * l2 * l1 + f * dl2_dx * l1 + f * l2 * dl1_dx);
  dgcon(0,2,2) = -(df_dx * l2 * l2 + f * dl2_dx * l2 + f * l2 * dl2_dx);
  dgcon(0,2,3) = -(df_dx * l2 * l3 + f * dl2_dx * l3 + f * l2 * dl3_dx);
  dgcon(0,3,0) = -(df_dx * l3 * l0 + f * dl3_dx * l0 + f * l3 * dl0_dx);
  dgcon(0,3,1) = -(df_dx * l3 * l1 + f * dl3_dx * l1 + f * l3 * dl1_dx);
  dgcon(0,3,2) = -(df_dx * l3 * l2 + f * dl3_dx * l2 + f * l3 * dl2_dx);
  dgcon(0,3,3) = -(df_dx * l3 * l3 + f * dl3_dx * l3 + f * l3 * dl3_dx);

  // Calculate metric component y-derivatives
  dgcon(1,0,0) = -(df_dy * l0 * l0 + f * dl0_dy * l0 + f * l0 * dl0_dy);
  dgcon(1,0,1) = -(df_dy * l0 * l1 + f * dl0_dy * l1 + f * l0 * dl1_dy);
  dgcon(1,0,2) = -(df_dy * l0 * l2 + f * dl0_dy * l2 + f * l0 * dl2_dy);
  dgcon(1,0,3) = -(df_dy * l0 * l3 + f * dl0_dy * l3 + f * l0 * dl3_dy);
  dgcon(1,1,0) = -(df_dy * l1 * l0 + f * dl1_dy * l0 + f * l1 * dl0_dy);
  dgcon(1,1,1) = -(df_dy * l1 * l1 + f * dl1_dy * l1 + f * l1 * dl1_dy);
  dgcon(1,1,2) = -(df_dy * l1 * l2 + f * dl1_dy * l2 + f * l1 * dl2_dy);
  dgcon(1,1,3) = -(df_dy * l1 * l3 + f * dl1_dy * l3 + f * l1 * dl3_dy);
  dgcon(1,2,0) = -(df_dy * l2 * l0 + f * dl2_dy * l0 + f * l2 * dl0_dy);
  dgcon(1,2,1) = -(df_dy * l2 * l1 + f * dl2_dy * l1 + f * l2 * dl1_dy);
  dgcon(1,2,2) = -(df_dy * l2 * l2 + f * dl2_dy * l2 + f * l2 * dl2_dy);
  dgcon(1,2,3) = -(df_dy * l2 * l3 + f * dl2_dy * l3 + f * l2 * dl3_dy);
  dgcon(1,3,0) = -(df_dy * l3 * l0 + f * dl3_dy * l0 + f * l3 * dl0_dy);
  dgcon(1,3,1) = -(df_dy * l3 * l1 + f * dl3_dy * l1 + f * l3 * dl1_dy);
  dgcon(1,3,2) = -(df_dy * l3 * l2 + f * dl3_dy * l2 + f * l3 * dl2_dy);
  dgcon(1,3,3) = -(df_dy * l3 * l3 + f * dl3_dy * l3 + f * l3 * dl3_dy);

  // Calculate metric component z-derivatives
  dgcon(2,0,0) = -(df_dz * l0 * l0 + f * dl0_dz * l0 + f * l0 * dl0_dz);
  dgcon(2,0,1) = -(df_dz * l0 * l1 + f * dl0_dz * l1 + f * l0 * dl1_dz);
  dgcon(2,0,2) = -(df_dz * l0 * l2 + f * dl0_dz * l2 + f * l0 * dl2_dz);
  dgcon(2,0,3) = -(df_dz * l0 * l3 + f * dl0_dz * l3 + f * l0 * dl3_dz);
  dgcon(2,1,0) = -(df_dz * l1 * l0 + f * dl1_dz * l0 + f * l1 * dl0_dz);
  dgcon(2,1,1) = -(df_dz * l1 * l1 + f * dl1_dz * l1 + f * l1 * dl1_dz);
  dgcon(2,1,2) = -(df_dz * l1 * l2 + f * dl1_dz * l2 + f * l1 * dl2_dz);
  dgcon(2,1,3) = -(df_dz * l1 * l3 + f * dl1_dz * l3 + f * l1 * dl3_dz);
  dgcon(2,2,0) = -(df_dz * l2 * l0 + f * dl2_dz * l0 + f * l2 * dl0_dz);
  dgcon(2,2,1) = -(df_dz * l2 * l1 + f * dl2_dz * l1 + f * l2 * dl1_dz);
  dgcon(2,2,2) = -(df_dz * l2 * l2 + f * dl2_dz * l2 + f * l2 * dl2_dz);
  dgcon(2,2,3) = -(df_dz * l2 * l3 + f * dl2_dz * l3 + f * l2 * dl3_dz);
  dgcon(2,3,0) = -(df_dz * l3 * l0 + f * dl3_dz * l0 + f * l3 * dl0_dz);
  dgcon(2,3,1) = -(df_dz * l3 * l1 + f * dl3_dz * l1 + f * l3 * dl1_dz);
  dgcon(2,3,2) = -(df_dz * l3 * l2 + f * dl3_dz * l2 + f * l3 * dl2_dz);
  dgcon(2,3,3) = -(df_dz * l3 * l3 + f * dl3_dz * l3 + f * l3 * dl3_dz);
  return;
}
