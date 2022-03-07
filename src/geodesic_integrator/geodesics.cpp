// Blacklight geodesic integrator - geodesic integration

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, ceil, isfinite, pow, sqrt
#include <sstream>    // ostringstream

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "geodesic_integrator.hpp"
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightWarning

//--------------------------------------------------------------------------------------------------

// Function for calculating ray positions and directions through space via Dormand-Prince
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes camera_pos[adaptive_level] and camera_dir[adaptive_level] have been set.
//   Initializes geodesic_num_steps[adaptive_level].
//   Allocates and initializes geodesic_pos, geodesic_dir, geodesic_len,
//       sample_flags[adaptive_level], and sample_num[adaptive_level].
//   Assumes x^0 is ignorable.
//   Integrates via the Dormand-Prince method (5th-order adaptive Runge-Kutta).
//     Method is RK5(4)7M of 1980 JCoAM 6 19.
//     4th-order interpolation follows 1986 MaCom 46 135, which gives a 4th-order estimate of the
//         midpoint and suggests combining this with values and derivatives at both endpoints to fit
//         a unique quartic over the step.
//     See Solving Ordinary Differential Equations I (Hairer, Norsett, Wanner) for coefficients that
//         accomplish the quartic fit without explicitly evaluating the midpoint.
//     All three references have different coefficients for the 4th-order step used in error
//         estimation; the coefficients here follow the original paper.
//     Interpolation is used to take steps small enough to satisfy user input ray_step, in that the
//         proper length of a step must be less than the product of ray_step with the radial
//         coordinate.
void GeodesicIntegrator::IntegrateGeodesicsDP()
{
  // Define coefficients
  double a_vals[7][6] = {};
  a_vals[1][0] = 1.0 / 5.0;
  a_vals[2][0] = 3.0 / 40.0;
  a_vals[2][1] = 9.0 / 40.0;
  a_vals[3][0] = 44.0 / 45.0;
  a_vals[3][1] = -56.0 / 15.0;
  a_vals[3][2] = 32.0 / 9.0;
  a_vals[4][0] = 19372.0 / 6561.0;
  a_vals[4][1] = -25360.0 / 2187.0;
  a_vals[4][2] = 64448.0 / 6561.0;
  a_vals[4][3] = -212.0 / 729.0;
  a_vals[5][0] = 9017.0 / 3168.0;
  a_vals[5][1] = -355.0 / 33.0;
  a_vals[5][2] = 46732.0 / 5247.0;
  a_vals[5][3] = 49.0 / 176.0;
  a_vals[5][4] = -5103.0 / 18656.0;
  a_vals[6][0] = 35.0 / 384.0;
  a_vals[6][2] = 500.0 / 1113.0;
  a_vals[6][3] = 125.0 / 192.0;
  a_vals[6][4] = -2187.0 / 6784.0;
  a_vals[6][5] = 11.0 / 84.0;
  double b_vals_5[7] =
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
  double b_vals_4[7] = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0,
      187.0 / 2100.0, 1.0 / 40.0};
  double b_vals_4m[7] = {6025192743.0 / 30085553152.0, 0.0, 51252292925.0 / 65400821598.0,
      -2691868925.0 / 45128329728.0, 187940372067.0 / 1594534317056.0,
      -1776094331.0 / 19743644256.0, 11237099.0 / 235043384.0};
  double d_vals[7] = {-12715105075.0 / 11282082432.0, 0.0, 87487479700.0 / 32700410799.0,
      -10690763975.0 / 1880347072.0, 701980252875.0 / 199316789632.0, -1453857185.0 / 822651844.0,
      69997945.0 / 29380423.0};

  // Define numerical parameters
  double err_power = 0.2;
  double ray_err_factor = 0.9;
  double ray_min_factor = 0.2;
  double ray_max_factor = 10.0;

  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  geodesic_pos.Allocate(num_pix, ray_max_steps, 4);
  geodesic_dir.Allocate(num_pix, ray_max_steps, 4);
  geodesic_len.Allocate(num_pix, ray_max_steps);
  geodesic_len.Zero();
  sample_flags[adaptive_level].Allocate(num_pix);
  sample_flags[adaptive_level].Zero();
  sample_num[adaptive_level].Allocate(num_pix);
  sample_num[adaptive_level].Zero();

  // Work in parallel
  int geodesic_num_steps_local = 0;
  int num_bad_geodesics = 0;
  #pragma omp parallel
  {
    // Allocate scratch arrays
    double gcon[4][4];
    double y_vals[9];
    double y_vals_temp[9];
    double y_vals_5[9];
    double y_vals_4[9];
    double y_vals_4m[8];
    double k_vals[7][9];
    double r_vals[4][8];

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      // Extract initial position
      y_vals[0] = camera_pos[adaptive_level](m,0);
      y_vals[1] = camera_pos[adaptive_level](m,1);
      y_vals[2] = camera_pos[adaptive_level](m,2);
      y_vals[3] = camera_pos[adaptive_level](m,3);

      // Extract initial momentum
      y_vals[4] = camera_dir[adaptive_level](m,0);
      y_vals[5] = camera_dir[adaptive_level](m,1);
      y_vals[6] = camera_dir[adaptive_level](m,2);
      y_vals[7] = camera_dir[adaptive_level](m,3);

      // Set initial proper distance
      y_vals[8] = 0.0;

      // Prepare to take steps
      for (int p = 0; p < 9; p++)
        y_vals_5[p] = y_vals[p];
      double r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);
      double h_new = -ray_step * r_new;
      int num_retry = 0;
      bool previous_fail = false;

      // Take steps
      for (int n = 0; n < ray_max_steps; )
      {
        // Check for too many retries
        if (num_retry > ray_max_retries)
        {
          sample_flags[adaptive_level](m) = true;
          break;
        }

        // Update step size
        double h = h_new;

        // Copy previous results
        if (not previous_fail and n > 0)
          for (int p = 0; p < 9; p++)
          {
            y_vals[p] = y_vals_5[p];
            k_vals[0][p] = k_vals[6][p];
          }
        if (not previous_fail and n == 0)
          GeodesicSubstepWithDistance(y_vals, k_vals[0]);
        double r = r_new;
        if (previous_fail)
          r = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);

        // Calculate substeps
        for (int substep = 1; substep < 7; substep++)
        {
          for (int p = 0; p < 9; p++)
            y_vals_temp[p] = y_vals[p];
          for (int q = 0; q < substep; q++)
            for (int p = 0; p < 9; p++)
              y_vals_temp[p] += a_vals[substep][q] * h * k_vals[q][p];
          GeodesicSubstepWithDistance(y_vals_temp, k_vals[substep]);
        }

        // Calculate values at end of full step
        for (int p = 0; p < 9; p++)
        {
          y_vals_5[p] = y_vals[p];
          y_vals_4[p] = y_vals[p];
        }
        for (int q = 0; q < 7; q++)
          for (int p = 0; p < 9; p++)
          {
            y_vals_5[p] += b_vals_5[q] * h * k_vals[q][p];
            y_vals_4[p] += b_vals_4[q] * h * k_vals[q][p];
          }
        r_new = RadialGeodesicCoordinate(y_vals_5[1], y_vals_5[2], y_vals_5[3]);

        // Estimate error
        double error = 0.0;
        for (int p = 0; p < 8; p++)
        {
          double y_abs = std::max(std::abs(y_vals[p]), std::abs(y_vals_5[p]));
          double error_scale = ray_tol_abs + ray_tol_rel * y_abs;
          double delta_y = std::abs(y_vals_5[p] - y_vals_4[p]);
          error = std::max(error, delta_y / error_scale);
        }

        // Decide if step is too far
        if (not (error <= 1.0))
        {
          double h_factor = ray_min_factor;
          if (std::isfinite(error))
          {
            double h_factor_ideal = ray_err_factor * std::pow(error, -err_power);
            h_factor = std::max(h_factor_ideal, ray_min_factor);
          }
          h_new = h * h_factor;
          num_retry += 1;
          previous_fail = true;
          continue;
        }
        else
        {
          double h_factor = ray_max_factor;
          if (error > 0.0)
          {
            h_factor = ray_err_factor * std::pow(error, -err_power);
            h_factor = std::max(h_factor, ray_min_factor);
            h_factor = std::min(h_factor, ray_max_factor);
          }
          if (previous_fail)
            h_factor = std::min(h_factor, 1.0);
          h_new = h * h_factor;
          num_retry = 0;
          previous_fail = false;
        }

        // Calculate values at middle of full step
        for (int p = 0; p < 8; p++)
          y_vals_4m[p] = y_vals[p];
        for (int q = 0; q < 7; q++)
          for (int p = 0; p < 8; p++)
            y_vals_4m[p] += b_vals_4m[q] * h * k_vals[q][p];

        // Subdivide full step
        double r_mid = RadialGeodesicCoordinate(y_vals_4m[1], y_vals_4m[2], y_vals_4m[3]);
        double delta_s_step = ray_step * r_mid;
        double delta_s_full = y_vals_5[8] - y_vals[8];
        int num_steps_ideal = static_cast<int>(std::ceil(delta_s_full / delta_s_step));
        delta_s_step = delta_s_full / num_steps_ideal;
        int num_steps_max = ray_max_steps - n;
        int num_steps = num_steps_ideal;
        if (num_steps > num_steps_max)
        {
          num_steps = num_steps_max;
          sample_flags[adaptive_level](m) = true;
        }

        // Calculate step midpoint if no subdivision necessary
        if (num_steps_ideal == 1)
        {
          geodesic_pos(m,n,0) = y_vals_4m[0];
          geodesic_pos(m,n,1) = y_vals_4m[1];
          geodesic_pos(m,n,2) = y_vals_4m[2];
          geodesic_pos(m,n,3) = y_vals_4m[3];
          geodesic_dir(m,n,0) = y_vals_4m[4];
          geodesic_dir(m,n,1) = y_vals_4m[5];
          geodesic_dir(m,n,2) = y_vals_4m[6];
          geodesic_dir(m,n,3) = y_vals_4m[7];
          geodesic_len(m,n) = h;
        }

        // Calculate interpolating coefficients for subdivisions
        if (num_steps_ideal > 1)
        {
          for (int p = 0; p < 8; p++)
          {
            r_vals[0][p] = y_vals_5[p] - y_vals[p];
            r_vals[1][p] = y_vals[p] - y_vals_5[p] + h * k_vals[0][p];
            r_vals[2][p] = 2.0 * (y_vals_5[p] - y_vals[p]) - h * (k_vals[0][p] + k_vals[6][p]);
            r_vals[3][p] = 0.0;
          }
          for (int q = 0; q < 7; q++)
            for (int p = 0; p < 8; p++)
              r_vals[3][p] += d_vals[q] * h * k_vals[q][p];
        }

        // Calculate subdivided steps
        if (num_steps_ideal > 1)
          for (int nn = 0; nn < num_steps; nn++)
          {
            double frac = (nn + 0.5) / num_steps_ideal;
            for (int p = 0; p < 8; p++)
              y_vals_temp[p] = y_vals[p] + frac * (r_vals[0][p] + (1.0 - frac) * (r_vals[1][p]
                  + frac * (r_vals[2][p] + (1.0 - frac) * r_vals[3][p])));
            geodesic_pos(m,n+nn,0) = y_vals_temp[0];
            geodesic_pos(m,n+nn,1) = y_vals_temp[1];
            geodesic_pos(m,n+nn,2) = y_vals_temp[2];
            geodesic_pos(m,n+nn,3) = y_vals_temp[3];
            geodesic_dir(m,n+nn,0) = y_vals_temp[4];
            geodesic_dir(m,n+nn,1) = y_vals_temp[5];
            geodesic_dir(m,n+nn,2) = y_vals_temp[6];
            geodesic_dir(m,n+nn,3) = y_vals_temp[7];
            geodesic_len(m,n+nn) = h / num_steps_ideal;
          }

        // Renormalize momentum
        ContravariantGeodesicMetric(y_vals_5[1], y_vals_5[2], y_vals_5[3], gcon);
        double temp_a = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_a += gcon[a][b] * y_vals_5[4+a] * y_vals_5[4+b];
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcon[0][a] * y_vals_5[4] * y_vals_5[4+a];
        double temp_c = gcon[0][0] * y_vals_5[4] * y_vals_5[4];
        double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
        double factor =
            temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
        for (int a = 1; a < 4; a++)
          y_vals_5[4+a] *= factor;

        // Check termination
        sample_num[adaptive_level](m) += num_steps;
        bool terminate_outer = r_new > camera_r and r_new > r;
        bool terminate_inner = r_new < r_terminate;
        if (terminate_outer or terminate_inner)
          break;
        bool last_step = n + num_steps >= ray_max_steps;
        if (last_step)
          sample_flags[adaptive_level](m) = true;

        // Prepare for next step
        n += num_steps;
      }
    }

    // Truncate geodesics at boundaries
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      int num_samples = sample_num[adaptive_level](m);
      if (num_samples > 1)
      {
        double r_new =
            RadialGeodesicCoordinate(geodesic_pos(m,0,1), geodesic_pos(m,0,2), geodesic_pos(m,0,3));
        for (int n = 1; n < num_samples; n++)
        {
          double r_old = r_new;
          r_new = RadialGeodesicCoordinate(geodesic_pos(m,n,1), geodesic_pos(m,n,2),
              geodesic_pos(m,n,3));
          bool terminate_outer = r_new > camera_r and r_new > r_old;
          bool terminate_inner = r_new < r_terminate;
          if (terminate_outer or terminate_inner)
          {
            sample_num[adaptive_level](m) = n;
            break;
          }
        }
      }
    }

    // Renormalize momenta
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
      for (int n = 0; n < sample_num[adaptive_level](m); n++)
      {
        ContravariantGeodesicMetric(geodesic_pos(m,n,1), geodesic_pos(m,n,2), geodesic_pos(m,n,3),
            gcon);
        double temp_a = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_a += gcon[a][b] * geodesic_dir(m,n,a) * geodesic_dir(m,n,b);
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcon[0][a] * geodesic_dir(m,n,0) * geodesic_dir(m,n,a);
        double temp_c = gcon[0][0] * geodesic_dir(m,n,0) * geodesic_dir(m,n,0);
        double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
        double factor =
            temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
        for (int a = 1; a < 4; a++)
          geodesic_dir(m,n,a) *= factor;
      }

    // Calculate maximum number of steps actually taken
    #pragma omp for schedule(static) reduction(max: geodesic_num_steps_local)
    for (int m = 0; m < num_pix; m++)
      geodesic_num_steps_local = std::max(geodesic_num_steps_local, sample_num[adaptive_level](m));

    // Calculate number of geodesics that do not terminate properly
    #pragma omp for schedule(static) reduction(+: num_bad_geodesics)
    for (int m = 0; m < num_pix; m++)
      if (sample_flags[adaptive_level](m))
        num_bad_geodesics++;
  }

  // Record number of steps taken
  geodesic_num_steps[adaptive_level] = geodesic_num_steps_local;

  // Report improperly terminated geodesics
  if (num_bad_geodesics > 0)
  {
    std::ostringstream message;
    message << num_bad_geodesics << " out of " << num_pix << " geodesics terminate unexpectedly.";
    BlacklightWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating ray positions and directions through space via 4th-order Runge-Kutta
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes camera_pos[adaptive_level] and camera_dir[adaptive_level] have been set.
//   Initializes geodesic_num_steps[adaptive_level].
//   Allocates and initializes geodesic_pos, geodesic_dir, geodesic_len,
//       sample_flags[adaptive_level], and sample_num[adaptive_level].
//   Assumes x^0 is ignorable.
//   Integrates via 4th-order Runge-Kutta with Butcher tableau
//        0  |
//       1/2 | 1/2
//       1/2 |  0  1/2
//        1  |  0   0   1
//            ----------------
//             1/6 1/3 1/3 1/6
//   Step size in affine parameter is taken to be the product of user input ray_step with the radial
//       displacement above the horizon.
void GeodesicIntegrator::IntegrateGeodesicsRK4()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  geodesic_pos.Allocate(num_pix, ray_max_steps, 4);
  geodesic_dir.Allocate(num_pix, ray_max_steps, 4);
  geodesic_len.Allocate(num_pix, ray_max_steps);
  geodesic_len.Zero();
  sample_flags[adaptive_level].Allocate(num_pix);
  sample_flags[adaptive_level].Zero();
  sample_num[adaptive_level].Allocate(num_pix);
  sample_num[adaptive_level].Zero();

  // Work in parallel
  int geodesic_num_steps_local = 0;
  int num_bad_geodesics = 0;
  #pragma omp parallel
  {
    // Allocate scratch arrays
    double gcon[4][4];
    double y_vals[8];
    double y_vals_substep[8];
    double y_vals_accumulate[8];
    double k_vals[8];

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      // Extract initial position
      y_vals[0] = camera_pos[adaptive_level](m,0);
      y_vals[1] = camera_pos[adaptive_level](m,1);
      y_vals[2] = camera_pos[adaptive_level](m,2);
      y_vals[3] = camera_pos[adaptive_level](m,3);
      double r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);

      // Extract initial momentum
      y_vals[4] = camera_dir[adaptive_level](m,0);
      y_vals[5] = camera_dir[adaptive_level](m,1);
      y_vals[6] = camera_dir[adaptive_level](m,2);
      y_vals[7] = camera_dir[adaptive_level](m,3);

      // Take steps
      for (int n = 0; n < ray_max_steps; n++)
      {
        // Calculate step size
        double r = r_new;
        double h = -ray_step * (r - r_horizon);

        // Calculate and accumulate first substep
        GeodesicSubstepWithoutDistance(y_vals, k_vals);
        for (int p = 0; p < 8; p++)
          y_vals_accumulate[p] = y_vals[p] + 1.0 / 6.0 * h * k_vals[p];

        // Calculate and accumulate second substep
        for (int p = 0; p < 8; p++)
          y_vals_substep[p] = y_vals[p] + 0.5 * h * k_vals[p];
        GeodesicSubstepWithoutDistance(y_vals_substep, k_vals);
        for (int p = 0; p < 8; p++)
          y_vals_accumulate[p] += 1.0 / 3.0 * h * k_vals[p];

        // Calculate and accumulate third substep
        for (int p = 0; p < 8; p++)
          y_vals_substep[p] = y_vals[p] + 0.5 * h * k_vals[p];
        GeodesicSubstepWithoutDistance(y_vals_substep, k_vals);
        for (int p = 0; p < 8; p++)
          y_vals_accumulate[p] += 1.0 / 3.0 * h * k_vals[p];

        // Calculate and accumulate fourth substep
        for (int p = 0; p < 8; p++)
          y_vals_substep[p] = y_vals[p] + h * k_vals[p];
        GeodesicSubstepWithoutDistance(y_vals_substep, k_vals);
        for (int p = 0; p < 8; p++)
          y_vals_accumulate[p] += 1.0 / 6.0 * h * k_vals[p];

        // Store midpoint
        for (int mu = 0; mu < 4; mu++)
        {
          geodesic_pos(m,n,mu) = 0.5 * (y_vals[mu] + y_vals_accumulate[mu]);
          geodesic_dir(m,n,mu) = 0.5 * (y_vals[4+mu] + y_vals_accumulate[4+mu]);
        }
        geodesic_len(m,n) = h;

        // Take step
        for (int p = 0; p < 8; p++)
          y_vals[p] = y_vals_accumulate[p];

        // Renormalize momentum
        ContravariantGeodesicMetric(y_vals[1], y_vals[2], y_vals[3], gcon);
        double temp_a = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_a += gcon[a][b] * y_vals[4+a] * y_vals[4+b];
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcon[0][a] * y_vals[4] * y_vals[4+a];
        double temp_c = gcon[0][0] * y_vals[4] * y_vals[4];
        double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
        double factor =
            temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
        for (int a = 1; a < 4; a++)
          y_vals[4+a] *= factor;

        // Check termination
        sample_num[adaptive_level](m)++;
        r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);
        bool terminate_outer = r_new > camera_r and r_new > r;
        bool terminate_inner = r_new < r_terminate;
        if (terminate_outer or terminate_inner)
          break;
        bool last_step = n + 1 >= ray_max_steps;
        if (last_step)
          sample_flags[adaptive_level](m) = true;
      }
    }

    // Truncate geodesics at boundaries
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      int num_samples = sample_num[adaptive_level](m);
      if (num_samples > 1)
      {
        double r_new =
            RadialGeodesicCoordinate(geodesic_pos(m,0,1), geodesic_pos(m,0,2), geodesic_pos(m,0,3));
        for (int n = 1; n < num_samples; n++)
        {
          double r_old = r_new;
          r_new = RadialGeodesicCoordinate(geodesic_pos(m,n,1), geodesic_pos(m,n,2),
              geodesic_pos(m,n,3));
          bool terminate_outer = r_new > camera_r and r_new > r_old;
          bool terminate_inner = r_new < r_terminate;
          if (terminate_outer or terminate_inner)
          {
            sample_num[adaptive_level](m) = n;
            break;
          }
        }
      }
    }

    // Renormalize momenta
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
      for (int n = 0; n < sample_num[adaptive_level](m); n++)
      {
        ContravariantGeodesicMetric(geodesic_pos(m,n,1), geodesic_pos(m,n,2), geodesic_pos(m,n,3),
            gcon);
        double temp_a = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_a += gcon[a][b] * geodesic_dir(m,n,a) * geodesic_dir(m,n,b);
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcon[0][a] * geodesic_dir(m,n,0) * geodesic_dir(m,n,a);
        double temp_c = gcon[0][0] * geodesic_dir(m,n,0) * geodesic_dir(m,n,0);
        double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
        double factor =
            temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
        for (int a = 1; a < 4; a++)
          geodesic_dir(m,n,a) *= factor;
      }

    // Calculate maximum number of steps actually taken
    #pragma omp for schedule(static) reduction(max: geodesic_num_steps_local)
    for (int m = 0; m < num_pix; m++)
      geodesic_num_steps_local = std::max(geodesic_num_steps_local, sample_num[adaptive_level](m));

    // Calculate number of geodesics that do not terminate properly
    #pragma omp for schedule(static) reduction(+: num_bad_geodesics)
    for (int m = 0; m < num_pix; m++)
      if (sample_flags[adaptive_level](m))
        num_bad_geodesics++;
  }

  // Record number of steps taken
  geodesic_num_steps[adaptive_level] = geodesic_num_steps_local;

  // Report improperly terminated geodesics
  if (num_bad_geodesics > 0)
  {
    std::ostringstream message;
    message << num_bad_geodesics << " out of " << num_pix << " geodesics terminate unexpectedly.";
    BlacklightWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating ray positions and directions through space via 2nd-order Runge-Kutta
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes camera_pos[adaptive_level] and camera_dir[adaptive_level] have been set.
//   Initializes geodesic_num_steps[adaptive_level].
//   Allocates and initializes geodesic_pos, geodesic_dir, geodesic_len,
//       sample_flags[adaptive_level], and sample_num[adaptive_level].
//   Assumes x^0 is ignorable.
//   Integrates via 2nd-order Runge-Kutta (Heun's method) with Butcher tableau
//       0 |
//       1 |  1
//          ----------------
//           1/2 1/2
//   Step size in affine parameter is taken to be the product of user input ray_step with the radial
//       displacement above the horizon.
void GeodesicIntegrator::IntegrateGeodesicsRK2()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  geodesic_pos.Allocate(num_pix, ray_max_steps, 4);
  geodesic_dir.Allocate(num_pix, ray_max_steps, 4);
  geodesic_len.Allocate(num_pix, ray_max_steps);
  geodesic_len.Zero();
  sample_flags[adaptive_level].Allocate(num_pix);
  sample_flags[adaptive_level].Zero();
  sample_num[adaptive_level].Allocate(num_pix);
  sample_num[adaptive_level].Zero();

  // Work in parallel
  int geodesic_num_steps_local = 0;
  int num_bad_geodesics = 0;
  #pragma omp parallel
  {
    // Allocate scratch arrays
    double gcon[4][4];
    double y_vals[8];
    double y_vals_substep[8];
    double k_vals[8];

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      // Extract initial position
      y_vals[0] = camera_pos[adaptive_level](m,0);
      y_vals[1] = camera_pos[adaptive_level](m,1);
      y_vals[2] = camera_pos[adaptive_level](m,2);
      y_vals[3] = camera_pos[adaptive_level](m,3);
      double r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);

      // Extract initial momentum
      y_vals[4] = camera_dir[adaptive_level](m,0);
      y_vals[5] = camera_dir[adaptive_level](m,1);
      y_vals[6] = camera_dir[adaptive_level](m,2);
      y_vals[7] = camera_dir[adaptive_level](m,3);

      // Take steps
      for (int n = 0; n < ray_max_steps; n++)
      {
        // Calculate step size
        double r = r_new;
        double h = -ray_step * (r - r_horizon);

        // Calculate and accumulate first substep
        GeodesicSubstepWithoutDistance(y_vals, k_vals);
        for (int p = 0; p < 8; p++)
          y_vals_substep[p] = y_vals[p] + h * k_vals[p];
        for (int p = 0; p < 8; p++)
          y_vals[p] += 1.0 / 2.0 * h * k_vals[p];

        // Store midpoint
        for (int mu = 0; mu < 4; mu++)
        {
          geodesic_pos(m,n,mu) = y_vals[mu];
          geodesic_dir(m,n,mu) = y_vals[4+mu];
        }
        geodesic_len(m,n) = h;

        // Calculate and accumulate second substep
        GeodesicSubstepWithoutDistance(y_vals_substep, k_vals);
        for (int p = 0; p < 8; p++)
          y_vals[p] += 1.0 / 2.0 * h * k_vals[p];

        // Renormalize momentum
        ContravariantGeodesicMetric(y_vals[1], y_vals[2], y_vals[3], gcon);
        double temp_a = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_a += gcon[a][b] * y_vals[4+a] * y_vals[4+b];
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcon[0][a] * y_vals[4] * y_vals[4+a];
        double temp_c = gcon[0][0] * y_vals[4] * y_vals[4];
        double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
        double factor =
            temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
        for (int a = 1; a < 4; a++)
          y_vals[4+a] *= factor;

        // Check termination
        sample_num[adaptive_level](m)++;
        r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);
        bool terminate_outer = r_new > camera_r and r_new > r;
        bool terminate_inner = r_new < r_terminate;
        if (terminate_outer or terminate_inner)
          break;
        bool last_step = n + 1 >= ray_max_steps;
        if (last_step)
          sample_flags[adaptive_level](m) = true;
      }
    }

    // Truncate geodesics at boundaries
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      int num_samples = sample_num[adaptive_level](m);
      if (num_samples > 1)
      {
        double r_new =
            RadialGeodesicCoordinate(geodesic_pos(m,0,1), geodesic_pos(m,0,2), geodesic_pos(m,0,3));
        for (int n = 1; n < num_samples; n++)
        {
          double r_old = r_new;
          r_new = RadialGeodesicCoordinate(geodesic_pos(m,n,1), geodesic_pos(m,n,2),
              geodesic_pos(m,n,3));
          bool terminate_outer = r_new > camera_r and r_new > r_old;
          bool terminate_inner = r_new < r_terminate;
          if (terminate_outer or terminate_inner)
          {
            sample_num[adaptive_level](m) = n;
            break;
          }
        }
      }
    }

    // Renormalize momenta
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
      for (int n = 0; n < sample_num[adaptive_level](m); n++)
      {
        ContravariantGeodesicMetric(geodesic_pos(m,n,1), geodesic_pos(m,n,2), geodesic_pos(m,n,3),
            gcon);
        double temp_a = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_a += gcon[a][b] * geodesic_dir(m,n,a) * geodesic_dir(m,n,b);
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcon[0][a] * geodesic_dir(m,n,0) * geodesic_dir(m,n,a);
        double temp_c = gcon[0][0] * geodesic_dir(m,n,0) * geodesic_dir(m,n,0);
        double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
        double factor =
            temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
        for (int a = 1; a < 4; a++)
          geodesic_dir(m,n,a) *= factor;
      }

    // Calculate maximum number of steps actually taken
    #pragma omp for schedule(static) reduction(max: geodesic_num_steps_local)
    for (int m = 0; m < num_pix; m++)
      geodesic_num_steps_local = std::max(geodesic_num_steps_local, sample_num[adaptive_level](m));

    // Calculate number of geodesics that do not terminate properly
    #pragma omp for schedule(static) reduction(+: num_bad_geodesics)
    for (int m = 0; m < num_pix; m++)
      if (sample_flags[adaptive_level](m))
        num_bad_geodesics++;
  }

  // Record number of steps taken
  geodesic_num_steps[adaptive_level] = geodesic_num_steps_local;

  // Report improperly terminated geodesics
  if (num_bad_geodesics > 0)
  {
    std::ostringstream message;
    message << num_bad_geodesics << " out of " << num_pix << " geodesics terminate unexpectedly.";
    BlacklightWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for reversing geodesics
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes geodesic_num_steps[adaptive_level], geodesic_pos, geodesic_dir, geodesic_len, and
//       sample_num[adaptive_level] have been set.
//   Allocates and initializes sample_pos[adaptive_level], sample_dir[adaptive_level], and
//       sample_len[adaptive_level], except reversed in the sampling dimension.
//   Deallocates geodesic_pos, geodesic_dir, and geodesic_len.
void GeodesicIntegrator::ReverseGeodesics()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  sample_pos[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level], 4);
  sample_dir[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level], 4);
  sample_len[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
  sample_len[adaptive_level].Zero();

  // Go through samples
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < num_pix; m++)
  {
    int num_steps = sample_num[adaptive_level](m);
    for (int n = 0; n < num_steps; n++)
    {
      // Skip terminated geodesics
      double len = geodesic_len(m,n);
      if (len == 0.0)
        break;

      // Set new arrays in reverse order
      sample_pos[adaptive_level](m,num_steps-1-n,0) = geodesic_pos(m,n,0);
      sample_pos[adaptive_level](m,num_steps-1-n,1) = geodesic_pos(m,n,1);
      sample_pos[adaptive_level](m,num_steps-1-n,2) = geodesic_pos(m,n,2);
      sample_pos[adaptive_level](m,num_steps-1-n,3) = geodesic_pos(m,n,3);
      sample_dir[adaptive_level](m,num_steps-1-n,0) = geodesic_dir(m,n,0);
      sample_dir[adaptive_level](m,num_steps-1-n,1) = geodesic_dir(m,n,1);
      sample_dir[adaptive_level](m,num_steps-1-n,2) = geodesic_dir(m,n,2);
      sample_dir[adaptive_level](m,num_steps-1-n,3) = geodesic_dir(m,n,3);
      sample_len[adaptive_level](m,num_steps-1-n) = -len;
    }
  }

  // Deallocate old arrays
  geodesic_pos.Deallocate();
  geodesic_dir.Deallocate();
  geodesic_len.Deallocate();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for taking single forward-Euler substep in time while computing proper distance
// Inputs:
//   y: dependent variables (positions, momenta, proper distance)
// Outputs:
//   k: derivatives with respect to independent variable (affine parameter)
// Notes:
//   Assumes k is allocated to be length 9.
//   Integrates the following equations:
//     d(x^mu) / d(lambda) = g^{mu nu} p_nu,
//     d(p_0) / d(lambda) = 0,
//     d(p_i) / d(lambda) = -1/2 * d(g^{mu nu}) / d(x^i) p_mu p_nu,
//     d(s) / d(lambda) = -(g_{i j} (g^{i mu} - g^{0 i} g^{0 mu} / g^{0 0}) p_mu
//         (g^{j nu} - g^{0 j} g^{0 nu} / g^{0 0}) p_nu)^(1/2).
//   Assumes x^0 is ignorable.
void GeodesicIntegrator::GeodesicSubstepWithDistance(double y[9], double k[9])
{
  double gcov[4][4];
  double gcon[4][4];
  double dgcon[3][4][4];
  CovariantGeodesicMetric(y[1], y[2], y[3], gcov);
  ContravariantGeodesicMetric(y[1], y[2], y[3], gcon);
  ContravariantGeodesicMetricDerivative(y[1], y[2], y[3], dgcon);
  for (int p = 0; p < 9; p++)
    k[p] = 0.0;
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      k[mu] += gcon[mu][nu] * y[4+nu];
  for (int a = 1; a < 4; a++)
    for (int mu = 0; mu < 4; mu++)
      for (int nu = 0; nu < 4; nu++)
        k[4+a] -= 0.5 * dgcon[a-1][mu][nu] * y[4+mu] * y[4+nu];
  double temp_a[4] = {};
  for (int a = 1; a < 4; a++)
    for (int mu = 0; mu < 4; mu++)
      temp_a[a] += (gcon[a][mu] - gcon[0][a] * gcon[0][mu] / gcon[0][0]) * y[4+mu];
  for (int a = 1; a < 4; a++)
    for (int b = 1; b < 4; b++)
      k[8] += gcov[a][b] * temp_a[a] * temp_a[b];
  k[8] = -std::sqrt(k[8]);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for taking single forward-Euler substep in time without computing proper distance
// Inputs:
//   y: dependent variables (positions, momenta)
// Outputs:
//   k: derivatives with respect to independent variable (affine parameter)
// Notes:
//   Assumes k is allocated to be length 8.
//   Integrates the following equations:
//     d(x^mu) / d(lambda) = g^{mu nu} p_nu,
//     d(p_0) / d(lambda) = 0,
//     d(p_i) / d(lambda) = -1/2 * d(g^{mu nu}) / d(x^i) p_mu p_nu,
//   Assumes x^0 is ignorable.
void GeodesicIntegrator::GeodesicSubstepWithoutDistance(double y[8], double k[8])
{
  double gcon[4][4];
  double dgcon[3][4][4];
  ContravariantGeodesicMetric(y[1], y[2], y[3], gcon);
  ContravariantGeodesicMetricDerivative(y[1], y[2], y[3], dgcon);
  for (int p = 0; p < 8; p++)
    k[p] = 0.0;
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      k[mu] += gcon[mu][nu] * y[4+nu];
  for (int a = 1; a < 4; a++)
    for (int mu = 0; mu < 4; mu++)
      for (int nu = 0; nu < 4; nu++)
        k[4+a] -= 0.5 * dgcon[a-1][mu][nu] * y[4+mu] * y[4+nu];
  return;
}
