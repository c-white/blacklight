// Blacklight ray tracer - geodesic integration

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, acos, atan, atan2, ceil, cos, hypot, isfinite, pow, sin, sqrt
#include <sstream>    // ostringstream

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "ray_tracer.hpp"
#include "../blacklight.hpp"        // math, enums
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightWarning

//--------------------------------------------------------------------------------------------------

// Function for setting up initial conditions for integrating geodesics
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes image_position and image_direction have been set except for the time components of
//       image_direction.
//   Initializes time components of image_direction.
//   Lowers all components of image_direction.
//   Quadratic solved as follows:
//     Outside ergosphere: unique positive root.
//     On ergosphere, assuming g_{0i} p^i < 0: unique root, which will be positive
//     Inside ergosphere, assuming g_{0i} p^i < 0: lesser positive root, which remains finite as
//         ergosphere is approached.
void RayTracer::InitializeGeodesics()
{
  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch array
    Array<double> gcov(4, 4);

    // Go through image pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        // Extract position
        double x[4];
        x[0] = image_position(m,l,0);
        x[1] = image_position(m,l,1);
        x[2] = image_position(m,l,2);
        x[3] = image_position(m,l,3);

        // Extract spatial components of momentum
        double p[4];
        p[1] = image_direction(m,l,1);
        p[2] = image_direction(m,l,2);
        p[3] = image_direction(m,l,3);

        // Calculate time component of momentum
        CovariantGeodesicMetric(x[1], x[2], x[3], gcov);
        double temp_a = gcov(0,0);
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcov(0,a) * p[a];
        double temp_c = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_c += gcov(a,b) * p[a] * p[b];
        double temp_d = std::sqrt(std::max(temp_b * temp_b - 4.0 * temp_a * temp_c, 0.0));
        p[0] = temp_a == 0 ? -temp_c / (2.0 * temp_b) : (temp_b < 0.0 ?
            2.0 * temp_c / (temp_d - temp_b) : -(temp_b + temp_d) / (2.0 * temp_a));

        // Lower momentum components
        for (int mu = 0; mu < 4; mu++)
        {
          image_direction(m,l,mu) = 0.0;
          for (int nu = 0; nu < 4; nu++)
            image_direction(m,l,mu) += gcov(mu,nu) * p[nu];
        }
      }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating ray positions and directions through space
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes image_position and image_direction have been set.
//   Initializes image_steps.
//   Allocates and initializes geodesic_pos, geodesic_dir, geodesic_len, sample_flags, and
//       sample_num.
//   Assumes x^0 is ignorable.
//   Integrates via the Dormand-Prince method (5th-order Runge-Kutta).
//     Method is RK5(4)7M of 1980 JCoAM 6 19.
//     4th-order interpolation follows 1986 MaCom 46 135, which gives a 4th-order estimate of the
//         midpoint and suggests combining this with values and derivatives at both endpoints to fit
//         a unique quartic over the step.
//     See Solving Ordinary Differential Equations I (Hairer, Norsett, Wanner) for coefficients that
//         accomplish the quartic fit without explicitly evaluating the midpoint.
//     All three references have different coefficients for the 4th-order step used in error
//         estimation; the coefficients here follow the original paper.
//     Interpolation is used to take steps small enough to satisfy user input ray_step.
void RayTracer::IntegrateGeodesics()
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

  // Define numerical parameter
  double err_power = 0.2;

  // Allocate arrays
  geodesic_pos.Allocate(image_resolution, image_resolution, ray_max_steps, 4);
  geodesic_dir.Allocate(image_resolution, image_resolution, ray_max_steps, 4);
  geodesic_len.Allocate(image_resolution, image_resolution, ray_max_steps);
  geodesic_len.Zero();
  sample_flags.Allocate(image_resolution, image_resolution);
  sample_flags.Zero();
  sample_num.Allocate(image_resolution, image_resolution);
  sample_num.Zero();

  // Work in parallel
  image_steps = 0;
  int num_bad_geodesics = 0;
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);
    Array<double> dgcon(3, 4, 4);
    double y_vals[9];
    double y_vals_temp[9];
    double y_vals_5[9];
    double y_vals_4[9];
    double y_vals_4m[8];
    double k_vals[7][9];
    double r_vals[4][8];

    // Go through image pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        // Extract initial position
        y_vals[0] = image_position(m,l,0);
        y_vals[1] = image_position(m,l,1);
        y_vals[2] = image_position(m,l,2);
        y_vals[3] = image_position(m,l,3);

        // Extract initial momentum
        y_vals[4] = image_direction(m,l,0);
        y_vals[5] = image_direction(m,l,1);
        y_vals[6] = image_direction(m,l,2);
        y_vals[7] = image_direction(m,l,3);

        // Set initial proper distance
        y_vals[8] = 0.0;

        // Prepare to take steps
        for (int p = 0; p < 9; p++)
          y_vals_5[p] = y_vals[p];
        double r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);
        double h_new = -r_new * ray_step;
        int num_retry = 0;
        bool previous_fail = false;

        // Take steps
        for (int n = 0; n < ray_max_steps; )
        {
          // Check for too many retries
          if (num_retry > ray_max_retries)
          {
            sample_flags(m,l) = true;
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
            for (int p = 0; p < 9; p++)
              GeodesicSubstep(y_vals, k_vals[0], gcov, gcon, dgcon);
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
            GeodesicSubstep(y_vals_temp, k_vals[substep], gcov, gcon, dgcon);
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
            sample_flags(m,l) = true;
          }

          // Calculate step midpoint if no subdivision necessary
          if (num_steps_ideal == 1)
          {
            geodesic_pos(m,l,n,0) = y_vals_4m[0];
            geodesic_pos(m,l,n,1) = y_vals_4m[1];
            geodesic_pos(m,l,n,2) = y_vals_4m[2];
            geodesic_pos(m,l,n,3) = y_vals_4m[3];
            geodesic_dir(m,l,n,0) = y_vals_4m[4];
            geodesic_dir(m,l,n,1) = y_vals_4m[5];
            geodesic_dir(m,l,n,2) = y_vals_4m[6];
            geodesic_dir(m,l,n,3) = y_vals_4m[7];
            geodesic_len(m,l,n) = h;
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
              geodesic_pos(m,l,n+nn,0) = y_vals_temp[0];
              geodesic_pos(m,l,n+nn,1) = y_vals_temp[1];
              geodesic_pos(m,l,n+nn,2) = y_vals_temp[2];
              geodesic_pos(m,l,n+nn,3) = y_vals_temp[3];
              geodesic_dir(m,l,n+nn,0) = y_vals_temp[4];
              geodesic_dir(m,l,n+nn,1) = y_vals_temp[5];
              geodesic_dir(m,l,n+nn,2) = y_vals_temp[6];
              geodesic_dir(m,l,n+nn,3) = y_vals_temp[7];
              geodesic_len(m,l,n+nn) = h / num_steps_ideal;
            }

          // Renormalize momentum
          ContravariantGeodesicMetric(y_vals_5[1], y_vals_5[2], y_vals_5[3], gcon);
          double temp_a = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              temp_a += gcon(a,b) * y_vals_5[4+a] * y_vals_5[4+b];
          double temp_b = 0.0;
          for (int a = 1; a < 4; a++)
            temp_b += 2.0 * gcon(0,a) * y_vals_5[4] * y_vals_5[4+a];
          double temp_c = gcon(0,0) * y_vals_5[4] * y_vals_5[4];
          double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
          double factor =
              temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
          for (int a = 1; a < 4; a++)
            y_vals_5[4+a] *= factor;

          // Check termination
          sample_num(m,l) += num_steps;
          bool terminate_outer = r_new > image_r and r_new > r;
          bool terminate_inner = r_new < r_terminate;
          if (terminate_outer or terminate_inner)
            break;
          bool last_step = n + num_steps >= ray_max_steps;
          if (last_step)
            sample_flags(m,l) = true;

          // Prepare for next step
          n += num_steps;
        }
      }

    // Truncate geodesics at boundaries
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        int num_samples = sample_num(m,l);
        if (num_samples > 1)
        {
          double r_new = RadialGeodesicCoordinate(geodesic_pos(m,l,0,1), geodesic_pos(m,l,0,2),
              geodesic_pos(m,l,0,3));
          for (int n = 1; n < num_samples; n++)
          {
            double r_old = r_new;
            r_new = RadialGeodesicCoordinate(geodesic_pos(m,l,n,1), geodesic_pos(m,l,n,2),
                geodesic_pos(m,l,n,3));
            bool terminate_outer = r_new > image_r and r_new > r_old;
            bool terminate_inner = r_new < r_terminate;
            if (terminate_outer or terminate_inner)
            {
              sample_num(m,l) = n;
              break;
            }
          }
        }
      }

    // Renormalize momenta
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
        for (int n = 0; n < sample_num(m,l); n++)
        {
          ContravariantGeodesicMetric(geodesic_pos(m,l,n,1), geodesic_pos(m,l,n,2),
              geodesic_pos(m,l,n,3), gcon);
          double temp_a = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              temp_a += gcon(a,b) * geodesic_dir(m,l,n,a) * geodesic_dir(m,l,n,b);
          double temp_b = 0.0;
          for (int a = 1; a < 4; a++)
            temp_b += 2.0 * gcon(0,a) * geodesic_dir(m,l,n,0) * geodesic_dir(m,l,n,a);
          double temp_c = gcon(0,0) * geodesic_dir(m,l,n,0) * geodesic_dir(m,l,n,0);
          double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
          double factor =
              temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
          for (int a = 1; a < 4; a++)
            geodesic_dir(m,l,n,a) *= factor;
        }

    // Calculate maximum number of steps actually taken
    #pragma omp for schedule(static) reduction(max:image_steps)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
        image_steps = std::max(image_steps, sample_num(m,l));

    // Calculate number of geodesics that do not terminate properly
    #pragma omp for schedule(static) reduction(+:num_bad_geodesics)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
        if (sample_flags(m,l))
          num_bad_geodesics++;
  }

  // Report improperly terminated geodesics
  if (num_bad_geodesics > 0)
  {
    std::ostringstream message;
    message << num_bad_geodesics << " out of " << image_resolution * image_resolution
        << " geodesics terminate unexpectedly.";
    BlacklightWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for transforming geodesics from integrating metric to simulation/formula metric
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes image_steps, geodesic_pos, geodesic_dir, geodesic_len, and sample_num have been set.
//   Allocates and initializes sample_pos, sample_dir, and sample_len, except transformed from
//       integrating metric to simulation/formula metric and reversed in the sampling dimension.
//   Deallocates geodesic_pos, geodesic_dir, and geodesic_len.
//   Assumes integrating and formula metrics are Cartesian Kerr-Schild.
//   Transformation of time components is trivial.
//   Transformation of geodesic length is trivial.
//   All transformations are trivial if both metrics are Cartesian Kerr-Schild.
void RayTracer::TransformGeodesics()
{
  // Allocate new arrays
  sample_pos.Allocate(image_resolution, image_resolution, image_steps, 4);
  sample_dir.Allocate(image_resolution, image_resolution, image_steps, 4);
  sample_len.Allocate(image_resolution, image_resolution, image_steps);
  sample_len.Zero();

  // Go through samples
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < image_resolution; m++)
    for (int l = 0; l < image_resolution; l++)
    {
      int num_steps = sample_num(m,l);
      for (int n = 0; n < num_steps; n++)
      {
        // Skip terminated geodesics
        double len = geodesic_len(m,l,n);
        if (len == 0.0)
          break;

        // Extract Cartesian position
        double t = geodesic_pos(m,l,n,0);
        double x = geodesic_pos(m,l,n,1);
        double y = geodesic_pos(m,l,n,2);
        double z = geodesic_pos(m,l,n,3);

        // Extract Cartesian direction
        double p_t = geodesic_dir(m,l,n,0);
        double p_x = geodesic_dir(m,l,n,1);
        double p_y = geodesic_dir(m,l,n,2);
        double p_z = geodesic_dir(m,l,n,3);

        // Prepare to find positions and momenta
        double x1, x2, x3;
        double p_1, p_2, p_3;

        // Account for model
        switch (model_type)
        {
          // Simulation output
          case ModelType::simulation:
          default:
          {
            // Account for simulation metric
            switch (simulation_coord)
            {
              // Spherical Kerr-Schild
              case Coordinates::sph_ks:
              default:
              {
                // Calculate spherical position
                double a2 = bh_a * bh_a;
                double rr2 = x * x + y * y + z * z;
                double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
                double r = std::sqrt(r2);
                double th = std::acos(z / r);
                double ph = std::atan2(y, x) - std::atan(bh_a / r);
                ph += ph < 0.0 ? 2.0 * math::pi : 0.0;
                ph -= ph >= 2.0 * math::pi ? 2.0 * math::pi : 0.0;
                double sth = std::sin(th);
                double cth = std::cos(th);
                double sph = std::sin(ph);
                double cph = std::cos(ph);

                // Calculate Jacobian of transformation
                double dx_dr = sth * cph;
                double dy_dr = sth * sph;
                double dz_dr = cth;
                double dx_dth = cth * (r * cph - bh_a * sph);
                double dy_dth = cth * (r * sph + bh_a * cph);
                double dz_dth = -r * sth;
                double dx_dph = sth * (-r * sph - bh_a * cph);
                double dy_dph = sth * (r * cph - bh_a * sph);
                double dz_dph = 0.0;

                // Calculate spherical direction
                double p_r = dx_dr * p_x + dy_dr * p_y + dz_dr * p_z;
                double p_th = dx_dth * p_x + dy_dth * p_y + dz_dth * p_z;
                double p_ph = dx_dph * p_x + dy_dph * p_y + dz_dph * p_z;

                // Assign position and direction
                x1 = r;
                x2 = th;
                x3 = ph;
                p_1 = p_r;
                p_2 = p_th;
                p_3 = p_ph;
                break;
              }

              // Cartesian Kerr-Schild
              case Coordinates::cart_ks:
              {
                x1 = x;
                x2 = y;
                x3 = z;
                p_1 = p_x;
                p_2 = p_y;
                p_3 = p_z;
                break;
              }
            }
            break;
          }

          // Formula
          case ModelType::formula:
          {
            x1 = x;
            x2 = y;
            x3 = z;
            p_1 = p_x;
            p_2 = p_y;
            p_3 = p_z;
            break;
          }
        }

        // Set new arrays in reverse order
        sample_pos(m,l,num_steps-1-n,0) = t;
        sample_pos(m,l,num_steps-1-n,1) = x1;
        sample_pos(m,l,num_steps-1-n,2) = x2;
        sample_pos(m,l,num_steps-1-n,3) = x3;
        sample_dir(m,l,num_steps-1-n,0) = p_t;
        sample_dir(m,l,num_steps-1-n,1) = p_1;
        sample_dir(m,l,num_steps-1-n,2) = p_2;
        sample_dir(m,l,num_steps-1-n,3) = p_3;
        sample_len(m,l,num_steps-1-n) = -len;
      }
    }

  // Deallocate old arrays
  geodesic_pos.Deallocate();
  geodesic_dir.Deallocate();
  geodesic_len.Deallocate();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for taking single forward-Euler substep in time
// Inputs:
//   y: dependent variables (positions, momenta, proper distance)
// Outputs:
//   k: derivatives with respect to independent variable (affine parameter)
//   gcov: components set
//   gcon: components set
//   dgcon: components set
// Notes:
//   Integrates following equations:
//     d(x^mu) / d(lambda) = g^{mu nu} p_nu
//     d(p_0) / d(lambda) = 0
//     d(p_i) / d(lambda) = -1/2 * d(g^{mu nu}) / d(x^i) p_mu p_nu
//     d(s) / d(lambda) = (g_{i j} g^{i mu} g^{j nu} p_mu p_nu)^(1/2)
//   Assumes x^0 is ignorable.
void RayTracer::GeodesicSubstep(double y[9], double k[9], Array<double> &gcov, Array<double> &gcon,
    Array<double> &dgcon)
{
  CovariantGeodesicMetric(y[1], y[2], y[3], gcov);
  ContravariantGeodesicMetric(y[1], y[2], y[3], gcon);
  ContravariantGeodesicMetricDerivative(y[1], y[2], y[3], dgcon);
  for (int p = 0; p < 9; p++)
    k[p] = 0.0;
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      k[mu] += gcon(mu,nu) * y[4+nu];
  for (int a = 1; a < 4; a++)
    for (int mu = 0; mu < 4; mu++)
      for (int nu = 0; nu < 4; nu++)
        k[4+a] -= 0.5 * dgcon(a-1,mu,nu) * y[4+mu] * y[4+nu];
  for (int a = 1; a < 4; a++)
    for (int b = 1; b < 4; b++)
      for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
          k[8] += gcov(a,b) * gcon(a,mu) * gcon(b,nu) * y[4+mu] * y[4+nu];
  k[8] = -std::sqrt(k[8]);
  return;
}
