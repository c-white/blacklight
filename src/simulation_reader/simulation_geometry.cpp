// Blacklight simulation reader - conversion functions for different coordinate systems

// C++ headers
#include <cmath>  // cos, exp, sin, sqrt

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "simulation_reader.hpp"
#include "../blacklight.hpp"      // Math
#include "../utils/array.hpp"     // Array

//--------------------------------------------------------------------------------------------------

// Function to convert coordinates from modified to standard spherical Kerr-Schild
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes x1f, x2f, x1v, and x2v are set in modified coordinates.
//   Allocates and sets x2v_alt to save modified coordinates.
//   Transforms x1f, x2f, x1v, and x2v.
//   Operates in serial, given all arrays are 1D and likely to be small.
void SimulationReader::ConvertCoordinates()
{
  // Extract parameters
  double h = metric_h;
  int n1 = x1v.n1;
  int n2 = x2v.n1;

  // Copy x^2 values
  x2v_alt.Allocate(n2);
  for (int j = 0; j < n2; j++)
    x2v_alt(j) = x2v(0,j);

  // Transform x^1 to r
  for (int i = 0; i <= n1; i++)
    x1f(0,i) = std::exp(x1f(0,i));
  for (int i = 0; i < n1; i++)
    x1v(0,i) = std::exp(x1v(0,i));

  // Transform x^2 to theta
  for (int j = 0; j <= n2; j++)
    x2f(0,j) = Math::pi * x2f(0,j) + (1.0 - h) / 2.0 * std::sin(2.0 * Math::pi * x2f(0,j));
  for (int j = 0; j < n2; j++)
    x2v(0,j) = Math::pi * x2v(0,j) + (1.0 - h) / 2.0 * std::sin(2.0 * Math::pi * x2v(0,j));
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to convert primitives from modified to standard spherical Kerr-Schild
// Inputs:
//   primitives: primitive array with modified vector components
// Outputs:
//   primitives: primitive array with standard vector components
// Notes:
//   Assumes x1v and x2v are set in standard coordinates.
//   Assumes x2v_alt is set in modified coordinates.
//   Assumes ind_uu1, ind_uu2, ind_uu3, ind_bb1, ind_bb2, and ind_bb3 are set.
void SimulationReader::ConvertPrimitives(Array<float> &primitives)
{
  // Extract parameters
  double a = simulation_a;
  double h = metric_h;
  int n1 = x1v.n1;
  int n2 = x2v.n1;
  int n3 = x3v.n1;

  // Go through cells in parallel
  #pragma omp parallel for schedule(static) collapse(2)
  for (int k = 0; k < n3; k++)
    for (int j = 0; j < n2; j++)
      for (int i = 0; i < n1; i++)
      {
        // Extract coordinates
        double r = x1v(0,i);
        double th = x2v(0,j);
        double sth = std::sin(th);
        double cth = std::cos(th);
        double x2 = x2v_alt(j);

        // Extract primitives
        double uu1 = primitives(ind_uu1,0,k,j,i);
        double uu2 = primitives(ind_uu2,0,k,j,i);
        double uu3 = primitives(ind_uu3,0,k,j,i);
        double bb1 = primitives(ind_bb1,0,k,j,i);
        double bb2 = primitives(ind_bb2,0,k,j,i);
        double bb3 = primitives(ind_bb3,0,k,j,i);

        // Calculate Jacobian of transformation
        double dr_dx1 = r;
        double dth_dx2 = Math::pi + (1.0 - h) * Math::pi * std::cos(2.0 * Math::pi * x2);

        // Calculate standard metric
        double sigma = r * r + a * a * cth * cth;
        double f = 2.0 * r / sigma;
        double g_tr = f;
        double g_tth = 0.0;
        double g_tph = -a * f * sth * sth;
        double g_rr = 1.0 + f;
        double g_rth = 0.0;
        double g_rph = -a * (1.0 + f) * sth * sth;
        double g_thth = sigma;
        double g_thph = 0.0;
        double g_phph = (r * r + a * a + a * a * f * sth * sth) * sth * sth;
        double gtt = -(1.0 + f);
        double gtr = f;
        double gtth = 0.0;
        double gtph = 0.0;
        double alpha = 1.0 / std::sqrt(-gtt);

        // Calculate modified metric
        double g_01 = dr_dx1 * g_tr;
        double g_02 = dth_dx2 * g_tth;
        double g_03 = g_tph;
        double g_11 = dr_dx1 * dr_dx1 * g_rr;
        double g_12 = dr_dx1 * dth_dx2 * g_rth;
        double g_13 = dr_dx1 * g_rph;
        double g_22 = dth_dx2 * dth_dx2 * g_thth;
        double g_23 = dth_dx2 * g_thph;
        double g_33 = g_phph;
        double g00 = gtt;
        double g01 = gtr / dr_dx1;
        double g02 = gtth / dth_dx2;
        double g03 = gtph;
        double alpha_mod = 1.0 / std::sqrt(-g00);

        // Transform velocity from modified normal frame to modified coordinate frame
        double uu0 = std::sqrt(1.0 + g_11 * uu1 * uu1 + 2.0 * g_12 * uu1 * uu2
            + 2.0 * g_13 * uu1 * uu3 + g_22 * uu2 * uu2 + 2.0 * g_23 * uu2 * uu3
            + g_33 * uu3 * uu3);
        double u0 = uu0 / alpha_mod;
        double u1 = uu1 - alpha_mod * g01 * uu0;
        double u2 = uu2 - alpha_mod * g02 * uu0;
        double u3 = uu3 - alpha_mod * g03 * uu0;
        double u_1 = g_01 * u0 + g_11 * u1 + g_12 * u2 + g_13 * u3;
        double u_2 = g_02 * u0 + g_12 * u1 + g_22 * u2 + g_23 * u3;
        double u_3 = g_03 * u0 + g_13 * u1 + g_23 * u2 + g_33 * u3;

        // Transform velocity from modified coordinate frame to standard coordinate frame
        double ut = u0;
        double ur = dr_dx1 * u1;
        double uth = dth_dx2 * u2;
        double uph = u3;

        // Transform velocity from standard coordinate frame to standard normal frame
        double uur = ur + alpha * alpha * gtr * ut;
        double uuth = uth + alpha * alpha * gtth * ut;
        double uuph = uph + alpha * alpha * gtph * ut;

        // Calculate magnetic 4-vector in modified coordinate frame
        double b0 = u_1 * bb1 + u_2 * bb2 + u_3 * bb3;
        double b1 = (bb1 + b0 * u1) / u0;
        double b2 = (bb2 + b0 * u2) / u0;
        double b3 = (bb3 + b0 * u3) / u0;

        // Transform magnetic 4-vector from modified coordinate frame to standard coordinate frame
        double bt = b0;
        double br = dr_dx1 * b1;
        double bth = dth_dx2 * b2;
        double bph = b3;

        // Calculate magnetic field in standard coordinate frame
        double bbr = br * ut - bt * ur;
        double bbth = bth * ut - bt * uth;
        double bbph = bph * ut - bt * uph;

        // Save results
        primitives(ind_uu1,0,k,j,i) = static_cast<float>(uur);
        primitives(ind_uu2,0,k,j,i) = static_cast<float>(uuth);
        primitives(ind_uu3,0,k,j,i) = static_cast<float>(uuph);
        primitives(ind_bb1,0,k,j,i) = static_cast<float>(bbr);
        primitives(ind_bb2,0,k,j,i) = static_cast<float>(bbth);
        primitives(ind_bb3,0,k,j,i) = static_cast<float>(bbph);
      }
  return;
}
