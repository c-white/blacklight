// Blacklight simulation reader - conversion functions for different coordinate systems

// C++ headers
#include <algorithm>  // min
#include <cmath>      // abs, cos, exp, log, pow, sin, sqrt

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "simulation_reader.hpp"
#include "../blacklight.hpp"        // Math
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function to convert coordinates from modified to standard spherical Kerr-Schild
// Inputs: (none)
// Outputs: (none)
// Notes:
//   If simulation coordinates are FMKS (MMKS):
//     Populates sks_map but leaves x[123][vf] unchanged (i.e., in native FMKS coordiantes).
//   Otherwise:
//     Assumes x1f, x2f, x1v, and x2v are set in modified coordinates.
//     Allocates and sets x2v_alt to save modified coordinates for regular MKS.
//     Transforms x1f, x2f, x1v, and x2v.
//     Operates in serial, given all arrays are 1D and likely to be small.
void SimulationReader::ConvertCoordinates()
{
  // Extract parameters
  double h = metric_h;
  int n1 = x1v.n1;
  int n2 = x2v.n1;

  // Handle FMKS case
  if (simulation_coord == Coordinates::fmks)
  {
    // Calculate radial limits
    double r_in = std::exp(x1f(0,0));
    double r_out = std::exp(x1f(0,n1));

    // Calculate map between x1/x2 and r/theta
    GenerateSKSMap(r_in, r_out);

    // Calculate grid bounds
    simulation_bounds.Allocate(6);
    double r_val = 0.0;
    double theta_val = 0.0;
    double phi_val = 0.0;
    GetSKSCoordinates(x1f(0,0), 0.0, 0.0, &r_val, &theta_val, &phi_val);
    simulation_bounds(0) = r_val;
    simulation_bounds(2) = theta_val;
    simulation_bounds(4) = phi_val;
    GetSKSCoordinates(x1f(0,n1), 1.0, 2.0 * Math::pi, &r_val, &theta_val, &phi_val);
    simulation_bounds(1) = r_val;
    simulation_bounds(3) = theta_val;
    simulation_bounds(5) = phi_val;
  }

  // Handle all other cases
  else
  {
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
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to convert primitive 3-vectors from modified to standard spherical Kerr-Schild
// Inputs:
//   primitives: primitive array with modified vector components
// Outputs:
//   primitives: primitive array with standard vector components
// Notes:
//   Assumes x1v and x2v are set in standard coordinates.
//   Assumes x2v_alt is set in modified coordinates.
//   Assumes ind_uu1, ind_uu2, ind_uu3, ind_bb1, ind_bb2, and ind_bb3 are set.
void SimulationReader::ConvertPrimitives3(Array<float> &primitives)
{
  // Extract parameters
  double a = simulation_a;
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
        double x1, x2;
        if (simulation_coord == Coordinates::sks)
        {
          x1 = std::log(r);
          x2 = x2v_alt(j);
        }
        else if (simulation_coord == Coordinates::fmks)
        {
          double phi;
          x1 = r;
          x2 = th;
          GetSKSCoordinates(x1, x2, x3v(0,k), &r, &th, &phi);
        }
        else
        {
          throw BlacklightException(
              "Attempting to translate primitives to CKS but could not process coordinates.");
        }
        double sth = std::sin(th);
        double cth = std::cos(th);

        // Extract primitives
        double uu1 = primitives(ind_uu1,0,k,j,i);
        double uu2 = primitives(ind_uu2,0,k,j,i);
        double uu3 = primitives(ind_uu3,0,k,j,i);
        double bb1 = primitives(ind_bb1,0,k,j,i);
        double bb2 = primitives(ind_bb2,0,k,j,i);
        double bb3 = primitives(ind_bb3,0,k,j,i);

        // Calculate Jacobian of transformation
        double dr_dx1, dth_dx1, dth_dx2;
        SetJacobianFactors(x1, x2, &dr_dx1, &dth_dx1, &dth_dx2);

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
        double g_01 = dr_dx1 * g_tr + dth_dx1 * g_tth;
        double g_02 = dth_dx2 * g_tth;
        double g_03 = g_tph;
        double g_11 =
            dr_dx1 * dr_dx1 * g_rr + 2.0 * dr_dx1 * dth_dx1 * g_rth + dth_dx1 * dth_dx1 * g_thth;
        double g_12 = dr_dx1 * dth_dx2 * g_rth + dth_dx1 * dth_dx2 * g_thth;
        double g_13 = dr_dx1 * g_rph + dth_dx1 * g_thph;
        double g_22 = dth_dx2 * dth_dx2 * g_thth;
        double g_23 = dth_dx2 * g_thph;
        double g_33 = g_phph;
        double g00 = gtt;
        double g01 = gtr / dr_dx1;
        double g02 = g_tth / dth_dx2 - dth_dx1 * g_tr / (dr_dx1 * dth_dx2);
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
        double uth = dth_dx1 * u1 + dth_dx2 * u2;
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
        double bth = dth_dx1 * b1 + dth_dx2 * b2;
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

//--------------------------------------------------------------------------------------------------

// Function to convert primitive 4-vectors from modified to standard spherical Kerr-Schild
// Inputs:
//   primitives: primitive array with modified vector components
// Outputs:
//   primitives: primitive array with standard vector components
// Notes:
//   Assumes x1v and x2v are set in standard coordinates.
//   Assumes x2v_alt is set in modified coordinates.
//   Assumes ind_u0, ind_uu1, ind_uu2, ind_uu3, ind_b0, ind_bb1, ind_bb2, and ind_bb3 are set.
void SimulationReader::ConvertPrimitives4(Array<float> &primitives)
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
        double cth = std::cos(th);
        double x2 = x2v_alt(j);

        // Extract primitives
        double u0 = primitives(ind_u0,0,k,j,i);
        double u1 = primitives(ind_uu1,0,k,j,i);
        double u2 = primitives(ind_uu2,0,k,j,i);
        double u3 = primitives(ind_uu3,0,k,j,i);
        double b0 = primitives(ind_b0,0,k,j,i);
        double b1 = primitives(ind_bb1,0,k,j,i);
        double b2 = primitives(ind_bb2,0,k,j,i);
        double b3 = primitives(ind_bb3,0,k,j,i);

        // Calculate Jacobian of transformation
        double dr_dx1 = r;
        double dth_dx2 = Math::pi + (1.0 - h) * Math::pi * std::cos(2.0 * Math::pi * x2);

        // Calculate standard metric
        double sigma = r * r + a * a * cth * cth;
        double f = 2.0 * r / sigma;
        double gtt = -(1.0 + f);
        double gtr = f;
        double gtth = 0.0;
        double gtph = 0.0;
        double alpha = 1.0 / std::sqrt(-gtt);

        // Transform velocity from modified coordinate frame to standard coordinate frame
        double ut = u0;
        double ur = dr_dx1 * u1;
        double uth = dth_dx2 * u2;
        double uph = u3;

        // Transform velocity from standard coordinate frame to standard normal frame
        double uur = ur + alpha * alpha * gtr * ut;
        double uuth = uth + alpha * alpha * gtth * ut;
        double uuph = uph + alpha * alpha * gtph * ut;

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

//--------------------------------------------------------------------------------------------------

// Function to generate the map between SKS and FMKS
// Inputs:
//   r_in: radial SKS coordinate for "inner edge" of coordinate map
//   r_out: radial SKS coordinate for "outer edge" of coordinate map
// Outputs: (none)
// Notes:
//   Assumes all metric parameters have been loaded.
//   Allocates and sets sks_map to save mapping.
//   Operates in serial since it only needs to run once and is 2D.
void SimulationReader::GenerateSKSMap(double r_in, double r_out)
{
  // TODO: these are not deallocated and will not work if the mesh coordinates change with time
  // Allocate map
  sks_map.Allocate(2, sks_map_n2, sks_map_n1);

  // Calculate spacing in SKS coordinates
  double dr = (r_out - r_in) / (sks_map_n1 - 1);
  double dtheta = Math::pi / (sks_map_n2 - 1);

  // Store parameters
  sks_map_r_in = r_in;
  sks_map_r_out = r_out;
  sks_map_dr = dr;
  sks_map_dtheta = dtheta;

  // Go through sample points in r
  for (int i = 0; i < sks_map_n1; ++i)
  {
    // Calculate radial coordinates
    double r = r_in + i * dr;
    double x1 = log(r);

    // Go through sample points in theta
    for (int j = 0; j < sks_map_n2; ++j)
    {
      // Calculate polar coordinates
      double theta = std::min(j * dtheta, Math::pi);
      double x2 = 0.5;

      // Iterate via bisection to find x^2
      if (theta > sks_map_tol and std::abs(Math::pi - theta) > sks_map_tol)
      {
        // Prepare initial bounds
        double x2_a = 0.0;
        double x2_b = 1.0;
        x2 = (x2_b + x2_a) / 2.0;
        double temp_r, temp_phi;
        double theta_a = 0.0;
        double theta_b = Math::pi;
        double theta_c = Math::pi / 2.0;
        GetSKSCoordinates(x1, x2_a, 0.0, &temp_r, &theta_a, &temp_phi);
        GetSKSCoordinates(x1, x2_b, 0.0, &temp_r, &theta_b, &temp_phi);

        // Perform iteration
        for (int n = 0; n < sks_map_max_iter; n++)
        {
          GetSKSCoordinates(x1, x2, 0.0, &temp_r, &theta_c, &temp_phi);
          if ((theta_c - theta) * (theta_b - theta) < 0.0)
          {
            theta_a = theta_c;
            x2_a = x2;
          }
          else
          {
            theta_b = theta_c;
            x2_b = x2;
          }
          x2 = (x2_a + x2_b) / 2.0;
          if (std::abs(theta - theta_c) < sks_map_tol)
            break;
        }
      }

      // Assign x^2 when at or beyond north pole
      else if (theta < sks_map_tol)
        x2 = 0.0;

      // Assign x^2 when at or beyond south pole
      else if (theta > Math::pi - sks_map_tol)
        x2 = 1.0;

      // Store mapping
      sks_map(0,j,i) = x1;
      sks_map(1,j,i) = x2;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to convert simulation coordinates to SKS
// Inputs:
//   x1, x2, x3: coordinate point in simulation coordinates
// Outputs:
//   *p_r, *p_theta, *p_phi: spherical Kerr-Schild coordinates for (x1, x2, x3)
void SimulationReader::GetSKSCoordinates(double x1, double x2, double x3, double *p_r,
    double *p_theta, double *p_phi)
{
  if (simulation_coord == Coordinates::fmks)
  {
    *p_r = std::exp(x1);
    double y = 2.0 * x2 - 1.0;
    double theta_g = Math::pi * x2 + (1.0 - metric_h) / 2.0 * std::sin(2.0 * Math::pi * x2);
    double theta_j = 0.5 * Math::pi + metric_derived_poly_norm * y
        * (1.0 + std::pow(y / metric_poly_xt, metric_poly_alpha) / (metric_poly_alpha + 1.0));
    *p_theta =
        theta_g + std::exp(metric_mks_smooth * (std::log(metric_r_in) - x1)) * (theta_j - theta_g);
    *p_phi = x3;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to calculate transformation factors from simulation coordinates to SKS
// Inputs:
//   x1, x2: simulation coordinates for where to compute the transformation
// Outputs:
//   *p_dr_dx1, *p_dth_dx1, *p_dth_dx2: Jacobian factors
void SimulationReader::SetJacobianFactors(double x1, double x2, double *p_dr_dx1, double *p_dth_dx1,
    double *p_dth_dx2)
{
  // Calculate dr/dx^1 in all cases
  *p_dr_dx1 = std::exp(x1);

  // Calculate theta elements of Jacobian in FMKS case
  if (simulation_coord == Coordinates::fmks)
  {
    double var_a = std::exp(metric_mks_smooth * (std::log(metric_r_in) - x1));
    double var_b = Math::pi * (0.5 - x2);
    double var_c = std::pow((2.0 * x2 - 1.0) / metric_poly_xt, metric_poly_alpha);
    double var_d = 1.0 + metric_poly_alpha;
    double var_e = metric_derived_poly_norm * (1.0 + var_c / var_d);
    double var_f = var_e * (2.0 * x2 - 1.0);
    double var_g = -0.5 * (1.0 - metric_h) * std::sin(2.0 * Math::pi * x2);
    *p_dth_dx1 = -metric_mks_smooth * var_a * (var_b + var_f + var_g);
    double var_h = Math::pi + (1.0 - metric_h) * Math::pi * std::cos(2.0 * Math::pi * x2);
    double var_i = -Math::pi + 2.0 * var_e;
    double var_j = 2.0 * metric_derived_poly_norm * metric_poly_alpha * var_c / var_d;
    double var_k = -(1.0 - metric_h) * Math::pi * std::cos(2.0 * Math::pi * x2);
    *p_dth_dx2 = var_h + var_a * (var_i + var_j + var_k);
  }

  // Calculate theta elements of Jacobian in MKS case
  else
  {
    *p_dth_dx1 = 0.0;
    *p_dth_dx2 = Math::pi + (1.0 - metric_h) * Math::pi * std::cos(2.0 * Math::pi * x2);
  }
  return;
}
