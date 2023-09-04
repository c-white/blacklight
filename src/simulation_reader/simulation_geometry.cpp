// Blacklight simulation reader - conversion functions for different coordinate systems

// C++ headers
#include <cmath>  // cos, exp, sin, sqrt

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
//     Populate sks_map but leave x[123][vf] unchanged, i.e., in native FMKS coordiantes.
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

  if (simulation_coord == Coordinates::fmks)
  {
    double r_in = std::exp(x1f(0, 0));
    double r_out = std::exp(x1f(0, n1));
    GenerateSKSMap(r_in, r_out, 2048, 2048);
    
    simulation_bounds.Allocate(6);

    double rval, thetaval, phival;
    GetSKSCoordinate(x1f(0, 0), 0.0, 0.0, rval, thetaval, phival);
    simulation_bounds(0) = rval;
    simulation_bounds(2) = thetaval;
    simulation_bounds(4) = phival;

    GetSKSCoordinate(x1f(0, n1), 1.0, 2.0*Math::pi, rval, thetaval, phival);
    simulation_bounds(1) = rval;
    simulation_bounds(3) = thetaval;
    simulation_bounds(5) = phival;

    // TODO: remove
    fprintf(stderr, "bounds: %g %g %g %g %g %g\n", simulation_bounds(0), simulation_bounds(1), simulation_bounds(2), simulation_bounds(3), simulation_bounds(4), simulation_bounds(5));
  }
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

// Function to return spatial SKS (r, theta, phi) coordinates given input coordiantes in
// native coordinate system.
// Inputs:
//   x1, x2, x3: coordinate point in simulation coordinates
// Outputs:
//   r, theta, phi: spherical Kerr-Schild coordinates for (x1, x2, x3)
void SimulationReader::GetSKSCoordinate(double x1, double x2, double x3, double &r, double &theta, double &phi)
{
  double h = metric_h;
  double poly_xt = metric_poly_xt;
  double poly_alpha = metric_poly_alpha;
  double mks_smooth = metric_mks_smooth;
  double rin = metric_rin;
  double poly_norm = metric_derived_poly_norm;

  if (simulation_coord == Coordinates::fmks) {
    r = exp(x1);
    phi = x3;
    double y = 2.0 * x2 - 1.0;
    double theta_G = Math::pi*x2 + ((1.0 - h) / 2.0) * std::sin(2.0 * Math::pi * x2);
    double theta_J = poly_norm * y * (1.0 + std::pow(y / poly_xt, poly_alpha) / (poly_alpha + 1.0));
    theta_J += 0.5 * Math::pi;
    theta = theta_G + std::exp(mks_smooth * (std::log(rin) - x1)) * (theta_J - theta_G);
  }
}

//--------------------------------------------------------------------------------------------------

// Function to generate the map between SKS (r, theta) coordinates and alternative
// spherically based coordinate systems like FMKS. Currently written for FMKS.
// Inputs: 
//   r_in: radial SKS coordinate for "inner edge" of coordinate map
//   r_out: radial SKS coordinate for "outer edge" of coordinate map
//   n1: number of SKS grid points (in radial direction) to sample over
//   n2: number of SKS grid points (in elevation/theta direction) to sample over
// Outputs: (none)
// Notes:
//   Assumes all metric parameters have been loaded.
//   Allocates and sets sks_map to save mapping.
//   Operates in serial since it only needs to run once and is 2D.
void SimulationReader::GenerateSKSMap(double r_in, double r_out, int n1, int n2)
{
  // TODO: these are not deallocated and will not work if the mesh coordinates change with time
  sks_map.Allocate(2, n2, n1);

  double dr = (r_out - r_in) / (n1 - 1);
  double dtheta = Math::pi / (n2 - 1);

  sks_map_rin = r_in;
  sks_map_rout = r_out;
  sks_map_dr = dr;
  sks_map_dtheta = dtheta;

  // TODO: this is a reasonable value, but maybe it should be set somewhere else
  double TOLERANCE = 1.e-8;

  for (int i = 0; i < n1; ++i)
  {
    double r = r_in + i * dr;
    double x1 = log(r);
    for (int j = 0; j < n2; ++j)
    {
      double theta = std::min(j * dtheta, Math::pi);
      double x2 = 0.5;

      if (theta > TOLERANCE && fabs(Math::pi - theta) > TOLERANCE)
      {
        // bisect down to correct value for x2
        double x2a = 0.0;
        double x2b = 1.0;
        x2 = (x2b + x2a) / 2.0;

        double tr, tphi;
        double theta_a = 0.0;
        double theta_b = Math::pi;
        double theta_c = Math::pi/2.0;
        GetSKSCoordinate(x1, x2a, 0.0, tr, theta_a, tphi);
        GetSKSCoordinate(x1, x2b, 0.0, tr, theta_b, tphi);

        for (int k = 0; k < 1000; ++k) {
          GetSKSCoordinate(x1, x2, 0.0, tr, theta_c, tphi);
          if ((theta_c - theta) * (theta_b - theta) < 0.0) {
            theta_a = theta_c;
            x2a = x2;
          } else {
            theta_b = theta_c;
            x2b = x2;
          }
          x2 = (x2a + x2b) / 2.0;
          if (fabs(theta - theta_c) < TOLERANCE) {
            break;
          }
        }
      } 
      
      else if (theta < TOLERANCE) 
      {
        x2 = 0.0;
      } 
      
      else if (fabs(Math::pi - theta) < TOLERANCE) 
      {
        x2 = 1.0;
      }

      // set coordinates in mesh
      sks_map(0, j, i) = x1;
      sks_map(1, j, i) = x2;
    }
  }

  return;
}

//--------------------------------------------------------------------------------------------------

// Function to set transformation factors from native (simulation) coordinates to spherical Kerr-Schild
// Inputs:
//   x1, x2 simulation coordinates for where to compute the transformation
// Outputs:
//   dr_dx1, dth_dx1, dth_dx2: Jacobian factors
void SimulationReader::SetJacobianFactors(double x1, double x2, double &dr_dx1, double &dth_dx1, double &dth_dx2)
{
  double h = metric_h;
  double poly_xt = metric_poly_xt;
  double poly_alpha = metric_poly_alpha;
  double mks_smooth = metric_mks_smooth;
  double rin = metric_rin;
  double poly_norm = metric_derived_poly_norm;

  // regular MKS
  dr_dx1 = std::exp(x1);
  dth_dx1 = 0.0;
  dth_dx2 = Math::pi + (1.0 - h) * Math::pi * std::cos(2.0 * Math::pi * x2);
  
  // FMKS
  if (simulation_coord == Coordinates::fmks)
  {
    dth_dx1 = - std::exp(mks_smooth * (std::log(rin) - x1)) * mks_smooth
            * (
            Math::pi / 2.0 -
            Math::pi * x2
                + poly_norm * (2.0 * x2 - 1.0)
                    * (1.0
                        + (std::pow((-1.0 + 2.0*x2) / poly_xt, poly_alpha))
                            / (1.0 + poly_alpha))
                - 1.0 / 2.0 * (1.0 - h) * std::sin(2.0 * Math::pi * x2));
    dth_dx2 = Math::pi + (1.0 - h) * Math::pi * std::cos(2.0 * Math::pi * x2)
            + std::exp(mks_smooth * (std::log(rin) - x1))
                * (-Math::pi
                    + 2.0 * poly_norm
                        * (1.0
                            + std::pow((2.0*x2 - 1.0) / poly_xt, poly_alpha)
                                / (poly_alpha + 1.0))
                    + (2.0 * poly_alpha * poly_norm * (2.0*x2 - 1.0)
                        * std::pow((2.0*x2 - 1.0) / poly_xt, poly_alpha - 1.0))
                        / ((1.0 + poly_alpha) * poly_xt)
                    - (1.0 - h) * Math::pi * std::cos(2.0 * Math::pi * x2));
  }
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
        double r, th, x1, x2;
        if (simulation_coord == Coordinates::sks)
        {
          r = x1v(0,i);
          th = x2v(0,j);
          x1 = std::log(r);
          x2 = x2v_alt(j);
        }
        else if (simulation_coord == Coordinates::fmks)
        {
          double phi;
          x1 = x1v(0,i);
          x2 = x2v(0,j);
          GetSKSCoordinate(x1, x2, x3v(0,k), r, th, phi);
        }
        else
        {
          throw BlacklightException("Attempting to translate primitives to CKS but could not process coordinates.");
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
        SetJacobianFactors(x1, x2, dr_dx1, dth_dx1, dth_dx2);

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
        double g_11 = dr_dx1*dr_dx1 * g_rr + 2.0 * dr_dx1 * dth_dx1 * g_rth + dth_dx1*dth_dx1 * g_thth;
        double g_12 = dr_dx1 * dth_dx2 * g_rth + dth_dx1 * dth_dx2 * g_thth;
        double g_13 = dr_dx1 * g_rph + dth_dx1 * g_thph;
        double g_22 = dth_dx2*dth_dx2 * g_thth;
        double g_23 = dth_dx2 * g_thph;
        double g_33 = g_phph;
        double g00 = gtt;
        double g01 = gtr / dr_dx1;
        double g02 = g_tth / dth_dx2 - (dth_dx1 * g_tr) / (dr_dx1 * dth_dx2);
        double g03 = gtph;
        double alpha_mod = 1.0 / std::sqrt(-g00);

        // Transform velocity from modified normal frame to modified coordinate frame
        // TODO check for missing terms .. does anything need g_00? 
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

        // TODO remove
        /*
        fprintf(stderr, "gcov %g %g %g %g %g %g %g %g %g\n", g_01, g_02, g_03, g_11, g_12, g_13, g_22, g_23, g_33);
        fprintf(stderr, "gcon %g %g %g %g\n", g00, g01, g02, g03);

        fprintf(stderr, "%d %d %d -> %g %g (%g %g)\n", i, j, k, r, th, x1, x2);
        fprintf(stderr, "%g %g %g\n", uu1, uu2, uu3);
        // prims are what we expect
        fprintf(stderr, "fmks %g %g %g %g %g %g %g\n", u0, u1, u2, u3, u_1, u_2, u_3);
        // ucon_fmks is what we expect
         */

        // Transform velocity from modified coordinate frame to standard coordinate frame
        double ut = u0;
        double ur = dr_dx1 * u1;
        double uth = dth_dx1 * u1 + dth_dx2 * u2;
        double uph = u3;

        // Transform velocity from standard coordinate frame to standard normal frame
        double uur = ur + alpha * alpha * gtr * ut;
        double uuth = uth + alpha * alpha * gtth * ut;
        double uuph = uph + alpha * alpha * gtph * ut;

        // TODO remove
        //fprintf(stderr, "ucon_ks %g %g %g %g\n", ut, ur, uth, uph);
        //fprintf(stderr, "uprim_ks %g %g %g\n", uur, uuth, uuph);
        // so this should be good ...
        //std::exit(4);

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

        // TODO remove
        /*
        fprintf(stderr, "fmks Bi bcon %g %g %g %g %g %g %g\n", bb1, bb2, bb3, b0, b1, b2, b3);
        fprintf(stderr, "bcon_ks %g %g %g %g\n", bt, br, bth, bph);
        fprintf(stderr, "bprim_ks %g %g %g\n", bbr, bbth, bbph);
        // these values also agree
        std::exit(4);
         */
        
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
  // TODO check that we are in mks rather than fmks?
  //      it seems unlikely to me that this would happen by accident but might be good to check.

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
