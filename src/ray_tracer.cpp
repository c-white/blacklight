// Blacklight ray tracer

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, acos, atan, atan2, cbrt, copysign, cos, cyl_bessel_k, exp, expm1,
                      // fmax, hypot, pow, sin, sqrt
#include <sstream>    // stringstream
#include <string>     // string

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "ray_tracer.hpp"
#include "array.hpp"        // Array
#include "blacklight.hpp"   // math, physics, enumerations
#include "exceptions.hpp"   // BlacklightWarning
#include "read_athena.hpp"  // AthenaReader
#include "read_input.hpp"   // InputReader

//--------------------------------------------------------------------------------------------------

// Ray tracer constructor
// Inputs:
//   input_reader: object containing input parameters read from input file
//   athena_reader: object containing raw data read from data file
RayTracer::RayTracer(const InputReader &input_reader, const AthenaReader &athena_reader)
{
  // Copy general input data
  model_type = input_reader.model_type;

  // Set parameters
  if (model_type == simulation)
  {
    bh_m = 1.0;
    bh_a = input_reader.simulation_a;
  }
  if (model_type == formula)
  {
    bh_m = 1.0;
    bh_a = input_reader.formula_spin;
  }

  // Copy simulation parameters
  if (model_type == simulation)
  {
    simulation_m_msun = input_reader.simulation_m_msun;
    simulation_rho_cgs = input_reader.simulation_rho_cgs;
    simulation_coord = input_reader.simulation_coord;
  }

  // Copy plasma parameters
  if (model_type == simulation)
  {
    plasma_mu = input_reader.plasma_mu;
    plasma_ne_ni = input_reader.plasma_ne_ni;
    plasma_rat_high = input_reader.plasma_rat_high;
    plasma_rat_low = input_reader.plasma_rat_low;
    plasma_sigma_max = input_reader.plasma_sigma_max;
  }

  // Copy formula parameters
  if (model_type == formula)
  {
    formula_mass = input_reader.formula_mass;
    formula_r0 = input_reader.formula_r0;
    formula_h = input_reader.formula_h;
    formula_l0 = input_reader.formula_l0;
    formula_q = input_reader.formula_q;
    formula_nup = input_reader.formula_nup;
    formula_cn0 = input_reader.formula_cn0;
    formula_alpha = input_reader.formula_alpha;
    formula_a = input_reader.formula_a;
    formula_beta = input_reader.formula_beta;
  }

  // Copy image parameters
  im_cam = input_reader.im_cam;
  im_r = input_reader.im_r;
  im_th = input_reader.im_th;
  im_ph = input_reader.im_ph;
  im_rot = input_reader.im_rot;
  im_width = input_reader.im_width;
  im_res = input_reader.im_res;
  im_freq = input_reader.im_freq;
  im_pole = input_reader.im_pole;

  // Copy ray-tracing parameters
  ray_step = input_reader.ray_step;
  ray_max_steps = input_reader.ray_max_steps;
  ray_sample_interp = input_reader.ray_sample_interp;
  ray_flat = input_reader.ray_flat;

  // Copy raw data scalars
  if (model_type == simulation)
  {
    x1_min = athena_reader.x1_min;
    x1_max = athena_reader.x1_max;
    x2_min = athena_reader.x2_min;
    x2_max = athena_reader.x2_max;
    x3_min = athena_reader.x3_min;
    x3_max = athena_reader.x3_max;
  }

  // Make shallow copies of raw data arrays
  if (model_type == simulation)
  {
    x1f = athena_reader.x1f;
    x2f = athena_reader.x2f;
    x3f = athena_reader.x3f;
    x1v = athena_reader.x1v;
    x2v = athena_reader.x2v;
    x3v = athena_reader.x3v;
    grid_rho = athena_reader.prim;
    grid_rho.Slice(5, athena_reader.ind_rho);
    grid_pgas = athena_reader.prim;
    grid_pgas.Slice(5, athena_reader.ind_pgas);
    grid_uu1 = athena_reader.prim;
    grid_uu1.Slice(5, athena_reader.ind_uu1);
    grid_uu2 = athena_reader.prim;
    grid_uu2.Slice(5, athena_reader.ind_uu2);
    grid_uu3 = athena_reader.prim;
    grid_uu3.Slice(5, athena_reader.ind_uu3);
    grid_bb1 = athena_reader.bb;
    grid_bb1.Slice(5, athena_reader.ind_bb1);
    grid_bb2 = athena_reader.bb;
    grid_bb2.Slice(5, athena_reader.ind_bb2);
    grid_bb3 = athena_reader.bb;
    grid_bb3.Slice(5, athena_reader.ind_bb3);
  }

  // Calculate inner photon orbit radius
  r_photon = 2.0 * bh_m * (1.0 + std::cos(2.0 / 3.0 * std::acos(-std::abs(bh_a) / bh_m)));
}

//--------------------------------------------------------------------------------------------------

// Ray tracer destructor
RayTracer::~RayTracer() {}

//--------------------------------------------------------------------------------------------------

// Top-level function for processing raw data into image
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes all data arrays have been set.
void RayTracer::MakeImage()
{
  InitializeCamera();
  InitializeGeodesics();
  IntegrateGeodesics();
  TransformGeodesics();
  if (model_type == simulation)
  {
    SampleSimulationAlongGeodesics();
    IntegrateSimulationRadiation();
  }
  if (model_type == formula)
    IntegrateFormulaRadiation();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for setting up camera pixels and initial ray directions
// Inputs: (none)
// Output: (none)
// Notes:
//   Allocates and initializes im_pos and im_dir except for time components of im_dir.
//   Neglects spacetime curvature at camera location.
//   Symbols:
//     n: unit outward normal
//     u: unit right vector
//     v: unit up vector
// TODO: calculate initial position, direction, and orientation more exactly
// TODO: account for flatness when doing above
void RayTracer::InitializeCamera()
{
  // Calculate trigonometric quantities
  double im_sth = std::sin(im_th);
  double im_cth = std::cos(im_th);
  double im_sph = std::sin(im_ph);
  double im_cph = std::cos(im_ph);
  double im_srot = std::sin(im_rot);
  double im_crot = std::cos(im_rot);

  // Calculate camera position
  double im_x = im_sth * (im_r * im_cph + bh_a * im_sph);
  double im_y = im_sth * (im_r * im_sph - bh_a * im_cph);
  double im_z = im_r * im_cth;

  // Calculate camera direction
  double im_nx = im_sth * im_cph;
  double im_ny = im_sth * im_sph;
  double im_nz = im_cth;

  // Calculate camera vertical orientation without rotation
  double up_x = 0.0;
  double up_y = 0.0;
  double up_z = 1.0;
  if (im_pole)
  {
    up_y = 1.0;
    up_z = 0.0;
  }
  double up_n = up_x * im_nx + up_y * im_ny + up_z * im_nz;
  double vx = up_x - up_n * im_nx;
  double vy = up_y - up_n * im_ny;
  double vz = up_z - up_n * im_nz;
  double v_norm = std::sqrt(vx * vx + vy * vy + vz * vz);
  vx /= v_norm;
  vy /= v_norm;
  vz /= v_norm;

  // Calculate camera horizontal orientation without rotation
  double ux = vy * im_nz - vz * im_ny;
  double uy = vz * im_nx - vx * im_nz;
  double uz = vx * im_ny - vy * im_nx;

  // Calculate camera orientation with rotation
  double im_ux = ux * im_crot - vx * im_srot;
  double im_uy = uy * im_crot - vy * im_srot;
  double im_uz = uz * im_crot - vz * im_srot;
  double im_vx = vx * im_crot + ux * im_srot;
  double im_vy = vy * im_crot + uy * im_srot;
  double im_vz = vz * im_crot + uz * im_srot;

  // Allocate arrays
  im_pos.Allocate(im_res, im_res, 4);
  im_dir.Allocate(im_res, im_res, 4);

  // Initialize arrays based on camera type
  switch (im_cam)
  {
    // Plane with parallel rays
    case plane:
    {
      #pragma omp parallel for schedule(static)
      for (int m = 0; m < im_res; m++)
        for (int l = 0; l < im_res; l++)
        {
          // Calculate position
          double u = (l - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
          double v = (m - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
          im_pos(m,l,0) = 0.0;
          im_pos(m,l,1) = im_x + u * im_ux + v * im_vx;
          im_pos(m,l,2) = im_y + u * im_uy + v * im_vy;
          im_pos(m,l,3) = im_z + u * im_uz + v * im_vz;

          // Set direction
          im_dir(m,l,1) = im_nx;
          im_dir(m,l,2) = im_ny;
          im_dir(m,l,3) = im_nz;
        }
      break;
    }

    // Point with converging rays
    case pinhole:
    default:
    {
      #pragma omp parallel for schedule(static)
      for (int m = 0; m < im_res; m++)
        for (int l = 0; l < im_res; l++)
        {
          // Set position
          im_pos(m,l,0) = 0.0;
          im_pos(m,l,1) = im_x;
          im_pos(m,l,2) = im_y;
          im_pos(m,l,3) = im_z;

          // Calculate direction
          double u = (l - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
          double v = (m - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
          double dir_x = im_x - (u * im_ux + v * im_vx);
          double dir_y = im_y - (u * im_uy + v * im_vy);
          double dir_z = im_z - (u * im_uz + v * im_vz);
          double dir_norm = std::hypot(dir_x, dir_y, dir_z);
          im_dir(m,l,1) = dir_x / dir_norm;
          im_dir(m,l,2) = dir_y / dir_norm;
          im_dir(m,l,3) = dir_z / dir_norm;
        }
      break;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for setting up initial conditions for integrating geodesics
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_pos and im_dir have been set except for time components of im_dir.
//   Initializes time components of im_dir.
//   Lowers all components of im_dir.
void RayTracer::InitializeGeodesics()
{
  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch array
    Array<double> gcov(4, 4);

    // Go through image pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Extract position
        double x[4];
        x[0] = im_pos(m,l,0);
        x[1] = im_pos(m,l,1);
        x[2] = im_pos(m,l,2);
        x[3] = im_pos(m,l,3);

        // Extract spatial components of momentum
        double p[4];
        p[1] = im_dir(m,l,1);
        p[2] = im_dir(m,l,2);
        p[3] = im_dir(m,l,3);

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
        double temp_q = -0.5 * (temp_b
            + std::copysign(1.0, temp_b) * std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c));
        p[0] = std::fmax(temp_q / temp_a, temp_c / temp_q);

        // Lower momentum components
        for (int mu = 0; mu < 4; mu++)
        {
          im_dir(m,l,mu) = 0.0;
          for (int nu = 0; nu < 4; nu++)
            im_dir(m,l,mu) += gcov(mu,nu) * p[nu];
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
//   Assumes im_pos and im_dir have been set.
//   Initializes im_steps.
//   Allocates and initializes geodesic_pos, geodesic_dir, geodesic_len, sample_flags, and
//       sample_num.
//   Assumes x^0 is ignorable.
//   Integrates via the midpoint method (2nd-order RK).
// TODO: calculate better step size
void RayTracer::IntegrateGeodesics()
{
  // Allocate arrays
  geodesic_pos.Allocate(im_res, im_res, ray_max_steps, 4);
  geodesic_dir.Allocate(im_res, im_res, ray_max_steps, 4);
  geodesic_len.Allocate(im_res, im_res, ray_max_steps);
  geodesic_len.Zero();
  sample_flags.Allocate(im_res, im_res);
  sample_flags.Zero();
  sample_num.Allocate(im_res, im_res);
  sample_num.Zero();

  // Work in parallel
  im_steps = 0;
  int num_bad_geodesics = 0;
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcon(4, 4);
    Array<double> dgcon(3, 4, 4);

    // Go through image pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Extract initial position
        double x[4];
        x[0] = im_pos(m,l,0);
        x[1] = im_pos(m,l,1);
        x[2] = im_pos(m,l,2);
        x[3] = im_pos(m,l,3);

        // Extract initial momentum
        double p[4];
        p[0] = im_dir(m,l,0);
        p[1] = im_dir(m,l,1);
        p[2] = im_dir(m,l,2);
        p[3] = im_dir(m,l,3);

        // Take steps
        for (int n = 0; n < ray_max_steps; n++)
        {
          // Calculate step size for going back to source
          double r = RadialGeodesicCoordinate(x[1], x[2], x[3]);
          double step = -ray_step * r;

          // Calculate position at half step, checking that step is worth taking
          ContravariantGeodesicMetric(x[1], x[2], x[3], gcon);
          double dx1[4] = {};
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              dx1[mu] += gcon(mu,nu) * p[nu];
          for (int mu = 0; mu < 4; mu++)
            geodesic_pos(m,l,n,mu) = x[mu] + step/2.0 * dx1[mu];
          double delta_r = RadialGeodesicCoordinate(geodesic_pos(m,l,n,1), geodesic_pos(m,l,n,2),
              geodesic_pos(m,l,n,3)) - r;
          if ((r > im_r and delta_r > 0.0) or r < r_photon)
            break;

          // Calculate momentum at half step
          ContravariantGeodesicMetricDerivative(x[1], x[2], x[3], dgcon);
          double dp1[4] = {};
          for (int a = 1; a <= 3; a++)
            for (int mu = 0; mu < 4; mu++)
              for (int nu = 0; nu < 4; nu++)
                dp1[a] -= 0.5 * dgcon(a-1,mu,nu) * p[mu] * p[nu];
          for (int mu = 0; mu < 4; mu++)
            geodesic_dir(m,l,n,mu) = p[mu] + step/2.0 * dp1[mu];

          // Renormalize momentum
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

          // Calculate position at full step
          double dx2[4] = {};
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              dx2[mu] += gcon(mu,nu) * geodesic_dir(m,l,n,nu);
          for (int mu = 0; mu < 4; mu++)
            x[mu] += step * dx2[mu];

          // Calculate momentum at full step
          ContravariantGeodesicMetricDerivative(geodesic_pos(m,l,n,1), geodesic_pos(m,l,n,2),
              geodesic_pos(m,l,n,3), dgcon);
          double dp2[4] = {};
          for (int a = 1; a <= 3; a++)
            for (int mu = 0; mu < 4; mu++)
              for (int nu = 0; nu < 4; nu++)
                dp2[a] -= 0.5 * dgcon(a-1,mu,nu) * geodesic_dir(m,l,n,mu) * geodesic_dir(m,l,n,nu);
          for (int mu = 0; mu < 4; mu++)
            p[mu] += step * dp2[mu];

          // Renormalize momentum
          ContravariantGeodesicMetric(x[1], x[2], x[3], gcon);
          temp_a = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              temp_a += gcon(a,b) * p[a] * p[b];
          temp_b = 0.0;
          for (int a = 1; a < 4; a++)
            temp_b += 2.0 * gcon(0,a) * p[0] * p[a];
          temp_c = gcon(0,0) * p[0] * p[0];
          temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
          factor =
              temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
          for (int a = 1; a < 4; a++)
            p[a] *= factor;

          // Store length of step
          geodesic_len(m,l,n) = -step;

          // Check for too many steps taken
          double r_new = RadialGeodesicCoordinate(x[1], x[2], x[3]);
          delta_r = r_new - r;
          if (n == ray_max_steps - 1 and not ((r > im_r and delta_r > 0.0) or r < r_photon))
            sample_flags(m,l) = true;
          sample_num(m,l)++;
        }
      }

    // Calculate maximum number of steps actually taken
    #pragma omp for schedule(static) reduction(max:im_steps)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        im_steps = std::max(im_steps, sample_num(m,l));

    // Calculate number of geodesics that do not terminate properly
    #pragma omp for schedule(static) reduction(+:num_bad_geodesics)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        if (sample_flags(m,l))
          num_bad_geodesics++;
  }

  // Report improperly terminated geodesics
  if (num_bad_geodesics > 0)
  {
    std::stringstream message;
    message << num_bad_geodesics << " out of " << im_res * im_res
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
//   Assumes im_steps, geodesic_pos, geodesic_dir, geodesic_len, and sample_num have been set.
//   Allocates and initializes sample_pos, sample_dir, and sample_len, except transformed from
//       integrating metric to simulation/formula metric and reversed in the sampling dimension.
//   Deallocates geodesic_pos, geodesic_dir, and geodesic_len.
//   Assumes integrating and formula metrics are Cartesian Kerr-Schild.
//   Transformation of time components is trivial.
//   Transformation of geodesic length is trivial.
//   All transformations trivial if both metrics are Cartesian Kerr-Schild.
void RayTracer::TransformGeodesics()
{
  // Allocate new arrays
  sample_pos.Allocate(im_res, im_res, im_steps, 4);
  sample_dir.Allocate(im_res, im_res, im_steps, 4);
  sample_len.Allocate(im_res, im_res, im_steps);
  sample_len.Zero();

  // Go through samples
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
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
          case simulation:
          {
            // Account for simulation metric
            switch (simulation_coord)
            {
              // Spherical Kerr-Schild
              case sph_ks:
              {
                // Calculate spherical position
                double a2 = bh_a * bh_a;
                double rr2 = x * x + y * y + z * z;
                double r2 =
                    0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
                double r = std::sqrt(r2);
                double th = std::acos(z / r);
                double ph = std::atan2(y, x) + std::atan(bh_a / r);
                ph += ph < 0.0 ? 2.0 * math::pi : 0.0;
                ph -= ph > 2.0 * math::pi ? 2.0 * math::pi : 0.0;
                double sth = std::sin(th);
                double cth = std::cos(th);
                double sph = std::sin(ph);
                double cph = std::cos(ph);

                // Calculate Jacobian of transformation
                double dx_dr = sth * cph;
                double dy_dr = sth * sph;
                double dz_dr = cth;
                double dx_dth = cth * (r * cph + bh_a * sph);
                double dy_dth = cth * (r * sph - bh_a * cph);
                double dz_dth = -r * sth;
                double dx_dph = sth * (-r * sph + bh_a * cph);
                double dy_dph = sth * (r * cph + bh_a * sph);
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
              case cart_ks:
              default:
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
          case formula:
          default:
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
        sample_len(m,l,num_steps-1-n) = len;
      }
    }

  // Deallocate old arrays
  geodesic_pos.Deallocate();
  geodesic_dir.Deallocate();
  geodesic_len.Deallocate();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for resampling simulation cell data onto rays.
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_steps, sample_flags, sample_num, sample_pos, and sample_len have been set.
//   Allocates and initializes sample_rho, sample_pgas, sample_uu1, sample_uu2, sample_uu3,
//       sample_bb1, sample_bb2, and sample_bb3.
//   If ray_sample_interp == false, uses primitives from cell containing geodesic sample point.
//   If ray_sample_interp == true, performs trilinear interpolation to geodesic sample point from
//       cell centers, using only data within the same block of cells (i.e. sometimes using
//       extrapolation).
//   TODO: interpolate across block boundaries
void RayTracer::SampleSimulationAlongGeodesics()
{
  // Allocate resampling arrays
  sample_rho.Allocate(im_res, im_res, im_steps);
  sample_pgas.Allocate(im_res, im_res, im_steps);
  sample_uu1.Allocate(im_res, im_res, im_steps);
  sample_uu2.Allocate(im_res, im_res, im_steps);
  sample_uu3.Allocate(im_res, im_res, im_steps);
  sample_bb1.Allocate(im_res, im_res, im_steps);
  sample_bb2.Allocate(im_res, im_res, im_steps);
  sample_bb3.Allocate(im_res, im_res, im_steps);

  // Work in parallel
  #pragma omp parallel
  {
    // Prepare bookkeeping
    int n_b = x1f.n2;
    int n_i = x1f.n1 - 1;
    int n_j = x2f.n1 - 1;
    int n_k = x3f.n1 - 1;
    int b = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    double x1_min_block = x1f(b,0);
    double x1_max_block = x1f(b,n_i);
    double x2_min_block = x2f(b,0);
    double x2_max_block = x2f(b,n_j);
    double x3_min_block = x3f(b,0);
    double x3_max_block = x3f(b,n_k);

    // Resample cell data onto geodesics
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Extract number of steps along this geodesic
        int num_steps = sample_num(m,l);

        // Set fallback values if geodesic poorly terminated
        if (sample_flags(m,l))
        {
          for (int n = 0; n < num_steps; n++)
          {
            sample_rho(m,l,n) = rho_fallback;
            sample_pgas(m,l,n) = pgas_fallback;
            sample_uu1(m,l,n) = uu1_fallback;
            sample_uu2(m,l,n) = uu2_fallback;
            sample_uu3(m,l,n) = uu3_fallback;
            sample_bb1(m,l,n) = bb1_fallback;
            sample_bb2(m,l,n) = bb2_fallback;
            sample_bb3(m,l,n) = bb3_fallback;
          }
          continue;
        }

        // Go along geodesic
        for (int n = 0; n < num_steps; n++)
        {
          // End if geodesic terminated
          if (sample_len(m,l,n) == 0.0)
            continue;

          // Extract coordinates
          double x1 = sample_pos(m,l,n,1);
          double x2 = sample_pos(m,l,n,2);
          double x3 = sample_pos(m,l,n,3);

          // Determine block
          if (x1 < x1_min_block or x1 > x1_max_block or x2 < x2_min_block or x2 > x2_max_block
              or x3 < x3_min_block or x3 > x3_max_block)
          {
            // Check if block contains position
            int b_new;
            double x1_min_temp = x1_min_block;
            double x1_max_temp = x1_max_block;
            double x2_min_temp = x2_min_block;
            double x2_max_temp = x2_max_block;
            double x3_min_temp = x3_min_block;
            double x3_max_temp = x3_max_block;
            for (b_new = 0; b_new < n_b; b_new++)
            {
              x1_min_temp = x1f(b_new,0);
              x1_max_temp = x1f(b_new,n_i);
              x2_min_temp = x2f(b_new,0);
              x2_max_temp = x2f(b_new,n_j);
              x3_min_temp = x3f(b_new,0);
              x3_max_temp = x3f(b_new,n_k);
              if (x1 >= x1_min_temp and x1 <= x1_max_temp and x2 >= x2_min_temp
                  and x2 <= x2_max_temp and x3 >= x3_min_temp and x3 <= x3_max_temp)
                break;
            }

            // Set fallback values if off grid
            if (b_new == n_b)
            {
              sample_rho(m,l,n) = rho_fallback;
              sample_pgas(m,l,n) = pgas_fallback;
              sample_uu1(m,l,n) = uu1_fallback;
              sample_uu2(m,l,n) = uu2_fallback;
              sample_uu3(m,l,n) = uu3_fallback;
              sample_bb1(m,l,n) = bb1_fallback;
              sample_bb2(m,l,n) = bb2_fallback;
              sample_bb3(m,l,n) = bb3_fallback;
              continue;
            }

            // Set newly found block as one to search
            b = b_new;
            x1_min_block = x1_min_temp;
            x1_max_block = x1_max_temp;
            x2_min_block = x2_min_temp;
            x2_max_block = x2_max_temp;
            x3_min_block = x3_min_temp;
            x3_max_block = x3_max_temp;
          }

          // Determine cell
          for (i = 0; i < n_i; i++)
            if (static_cast<double>(x1f(b,i+1)) >= x1)
              break;
          for (j = 0; j < n_j; j++)
            if (static_cast<double>(x2f(b,j+1)) >= x2)
              break;
          for (k = 0; k < n_k; k++)
            if (static_cast<double>(x3f(b,k+1)) >= x3)
              break;

          // Resample values without interpolation
          if (not ray_sample_interp)
          {
            sample_rho(m,l,n) = grid_rho(b,k,j,i);
            sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
            sample_uu1(m,l,n) = grid_uu1(b,k,j,i);
            sample_uu2(m,l,n) = grid_uu2(b,k,j,i);
            sample_uu3(m,l,n) = grid_uu3(b,k,j,i);
            sample_bb1(m,l,n) = grid_bb1(b,k,j,i);
            sample_bb2(m,l,n) = grid_bb2(b,k,j,i);
            sample_bb3(m,l,n) = grid_bb3(b,k,j,i);
          }

          // Resample values with interpolation
          if (ray_sample_interp)
          {

            // Calculate interpolation/extrapolation indices and coefficients
            int i_m = i == 0 or (i != n_i - 1 and x1 >= static_cast<double>(x1v(b,i))) ? i : i - 1;
            int i_p = i_m + 1;
            int j_m = j == 0 or (j != n_j - 1 and x2 >= static_cast<double>(x2v(b,j))) ? j : j - 1;
            int j_p = j_m + 1;
            int k_m = k == 0 or (k != n_k - 1 and x3 >= static_cast<double>(x3v(b,k))) ? k : k - 1;
            int k_p = k_m + 1;
            double f_i = (x1 - static_cast<double>(x1v(b,i_m)))
                / (static_cast<double>(x1v(b,i_p)) - static_cast<double>(x1v(b,i_m)));
            double f_j = (x2 - static_cast<double>(x2v(b,j_m)))
                / (static_cast<double>(x2v(b,j_p)) - static_cast<double>(x2v(b,j_m)));
            double f_k = (x3 - static_cast<double>(x3v(b,k_m)))
                / (static_cast<double>(x3v(b,k_p)) - static_cast<double>(x3v(b,k_m)));

            // Interpolate density
            double val_mmm = static_cast<double>(grid_rho(b,k_m,j_m,i_m));
            double val_mmp = static_cast<double>(grid_rho(b,k_m,j_m,i_p));
            double val_mpm = static_cast<double>(grid_rho(b,k_m,j_p,i_m));
            double val_mpp = static_cast<double>(grid_rho(b,k_m,j_p,i_p));
            double val_pmm = static_cast<double>(grid_rho(b,k_p,j_m,i_m));
            double val_pmp = static_cast<double>(grid_rho(b,k_p,j_m,i_p));
            double val_ppm = static_cast<double>(grid_rho(b,k_p,j_p,i_m));
            double val_ppp = static_cast<double>(grid_rho(b,k_p,j_p,i_p));
            sample_rho(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);
            if (sample_rho(m,l,n) <= 0.0f)
              sample_rho(m,l,n) = grid_rho(b,k,j,i);

            // Interpolate gas pressure
            val_mmm = static_cast<double>(grid_pgas(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_pgas(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_pgas(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_pgas(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_pgas(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_pgas(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_pgas(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_pgas(b,k_p,j_p,i_p));
            sample_pgas(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);
            if (sample_pgas(m,l,n) <= 0.0f)
              sample_pgas(m,l,n) = grid_pgas(b,k,j,i);

            // Interpolate x1-velocity
            val_mmm = static_cast<double>(grid_uu1(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_uu1(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_uu1(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_uu1(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_uu1(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_uu1(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_uu1(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_uu1(b,k_p,j_p,i_p));
            sample_uu1(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x2-velocity
            val_mmm = static_cast<double>(grid_uu2(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_uu2(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_uu2(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_uu2(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_uu2(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_uu2(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_uu2(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_uu2(b,k_p,j_p,i_p));
            sample_uu2(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x3-velocity
            val_mmm = static_cast<double>(grid_uu3(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_uu3(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_uu3(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_uu3(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_uu3(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_uu3(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_uu3(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_uu3(b,k_p,j_p,i_p));
            sample_uu3(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x1-field
            val_mmm = static_cast<double>(grid_bb1(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_bb1(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_bb1(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_bb1(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_bb1(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_bb1(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_bb1(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_bb1(b,k_p,j_p,i_p));
            sample_bb1(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x2-field
            val_mmm = static_cast<double>(grid_bb2(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_bb2(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_bb2(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_bb2(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_bb2(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_bb2(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_bb2(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_bb2(b,k_p,j_p,i_p));
            sample_bb2(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x3-field
            val_mmm = static_cast<double>(grid_bb3(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_bb3(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_bb3(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_bb3(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_bb3(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_bb3(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_bb3(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_bb3(b,k_p,j_p,i_p));
            sample_bb3(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);
          }
        }
      }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on simulation data
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_dir, sample_len, sample_rho, sample_pgas, sample_uu1, sample_uu2,
//       sample_uu3, sample_bb1, sample_bb2, and sample_bb3 have been set.
//   Allocates and initializes image.
//   Assumes x^0 is ignorable.
//   References symphony paper 2016 ApJ 822 34 (S).
void RayTracer::IntegrateSimulationRadiation()
{
  // Allocate image array
  image.Allocate(im_res, im_res);
  image.Zero();

  // Calculate units
  double x_unit = physics::gg_msun * simulation_m_msun / (physics::c * physics::c);
  double e_unit = simulation_rho_cgs * physics::c * physics::c;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Go through samples
        int num_steps = sample_num(m,l);
        for (int n = 0; n < num_steps; n++)
        {
          // Check that this sample contributes
          double delta_lambda = sample_len(m,l,n);
          if (delta_lambda == 0.0)
            continue;

          // Extract geodesic position and momentum
          double x1 = sample_pos(m,l,n,1);
          double x2 = sample_pos(m,l,n,2);
          double x3 = sample_pos(m,l,n,3);
          double p_0 = sample_dir(m,l,n,0);
          double p_1 = sample_dir(m,l,n,1);
          double p_2 = sample_dir(m,l,n,2);
          double p_3 = sample_dir(m,l,n,3);

          // Extract model variables
          double rho = sample_rho(m,l,n);
          double pgas = sample_pgas(m,l,n);
          double uu1 = sample_uu1(m,l,n);
          double uu2 = sample_uu2(m,l,n);
          double uu3 = sample_uu3(m,l,n);
          double bb1 = sample_bb1(m,l,n);
          double bb2 = sample_bb2(m,l,n);
          double bb3 = sample_bb3(m,l,n);

          // Skip contribution if magnetic field vanishes
          if (bb1 == 0.0 and bb2 == 0.0 and bb3 == 0.0)
            continue;

          // Calculate metric
          CovariantCoordinateMetric(x1, x2, x3, gcov);
          ContravariantCoordinateMetric(x1, x2, x3, gcon);

          // Calculate frequency in CGS units
          double p0 = gcon(0,0) * p_0 + gcon(0,1) * p_1 + gcon(0,2) * p_2 + gcon(0,3) * p_3;
          double nu_cgs = -p0 / p_0 * im_freq;

          // Calculate 4-velocity
          double temp = gcov(1,1) * uu1 * uu1 + 2.0 * gcov(1,2) * uu1 * uu2
              + 2.0 * gcov(1,3) * uu1 * uu3 + gcov(2,2) * uu2 * uu2 + 2.0 * gcov(2,3) * uu2 * uu3
              + gcov(3,3) * uu3 * uu3;
          double gamma = std::sqrt(1.0 + temp);
          double alpha = 1.0 / std::sqrt(-gcon(0,0));
          double beta1 = -gcon(0,1) / gcon(0,0);
          double beta2 = -gcon(0,2) / gcon(0,0);
          double beta3 = -gcon(0,3) / gcon(0,0);
          double u0 = gamma / alpha;
          double u1 = uu1 - beta1 * gamma / alpha;
          double u2 = uu2 - beta2 * gamma / alpha;
          double u3 = uu3 - beta3 * gamma / alpha;
          double u_0 = gcov(0,0) * u0 + gcov(0,1) * u1 + gcov(0,2) * u2 + gcov(0,3) * u3;
          double u_1 = gcov(1,0) * u0 + gcov(1,1) * u1 + gcov(1,2) * u2 + gcov(1,3) * u3;
          double u_2 = gcov(2,0) * u0 + gcov(2,1) * u1 + gcov(2,2) * u2 + gcov(2,3) * u3;
          double u_3 = gcov(3,0) * u0 + gcov(3,1) * u1 + gcov(3,2) * u2 + gcov(3,3) * u3;

          // Calculate magnetic field strength
          double b0 = u_1 * bb1 + u_2 * bb2 + u_3 * bb3;
          double b1 = (bb1 + b0 * u1) / u0;
          double b2 = (bb2 + b0 * u2) / u0;
          double b3 = (bb3 + b0 * u3) / u0;
          double b_0 = gcov(0,0) * b0 + gcov(0,1) * b1 + gcov(0,2) * b2 + gcov(0,3) * b3;
          double b_1 = gcov(1,0) * b0 + gcov(1,1) * b1 + gcov(1,2) * b2 + gcov(1,3) * b3;
          double b_2 = gcov(2,0) * b0 + gcov(2,1) * b1 + gcov(2,2) * b2 + gcov(2,3) * b3;
          double b_3 = gcov(3,0) * b0 + gcov(3,1) * b1 + gcov(3,2) * b2 + gcov(3,3) * b3;
          double b_sq = b_0 * b0 + b_1 * b1 + b_2 * b2 + b_3 * b3;

          // Skip contribution if considered to be vacuum
          if (plasma_sigma_max >= 0.0 and b_sq / rho > plasma_sigma_max)
            continue;

          // Calculate fluid-frame photon momentum
          double p_1f = p_1 - u_1 / u_0 * p_0;
          double p_2f = p_2 - u_2 / u_0 * p_0;
          double p_3f = p_3 - u_3 / u_0 * p_0;
          double g1f1f = gcon(1,1) + u1 * u1;
          double g1f2f = gcon(1,2) + u1 * u2;
          double g1f3f = gcon(1,3) + u1 * u3;
          double g2f2f = gcon(2,2) + u2 * u2;
          double g2f3f = gcon(2,3) + u2 * u3;
          double g3f3f = gcon(3,3) + u3 * u3;
          double p1f = g1f1f * p_1f + g1f2f * p_2f + g1f3f * p_3f;
          double p2f = g1f2f * p_1f + g2f2f * p_2f + g2f3f * p_3f;
          double p3f = g1f3f * p_1f + g2f3f * p_2f + g3f3f * p_3f;

          // Calculate fluid-frame magnetic field
          double bb1f = b1;
          double bb2f = b2;
          double bb3f = b3;
          double g_1f1f =
              gcov(1,1) - 2.0 * u_1 / u_0 * gcov(0,1) + u_1 * u_1 / (u_0 * u_0) * gcov(0,0);
          double g_1f2f = gcov(1,2) - u_1 / u_0 * gcov(0,2) - u_2 / u_0 * gcov(0,1)
              + u_1 * u_2 / (u_0 * u_0) * gcov(0,0);
          double g_1f3f = gcov(1,3) - u_1 / u_0 * gcov(0,3) - u_3 / u_0 * gcov(0,1)
              + u_1 * u_3 / (u_0 * u_0) * gcov(0,0);
          double g_2f2f =
              gcov(2,2) - 2.0 * u_2 / u_0 * gcov(0,2) + u_2 * u_2 / (u_0 * u_0) * gcov(0,0);
          double g_2f3f = gcov(2,3) - u_2 / u_0 * gcov(0,3) - u_3 / u_0 * gcov(0,2)
              + u_2 * u_3 / (u_0 * u_0) * gcov(0,0);
          double g_3f3f =
              gcov(3,3) - 2.0 * u_3 / u_0 * gcov(0,3) + u_3 * u_3 / (u_0 * u_0) * gcov(0,0);
          double bb_1f = g_1f1f * bb1f + g_1f2f * bb2f + g_1f3f * bb3f;
          double bb_2f = g_1f2f * bb1f + g_2f2f * bb2f + g_2f3f * bb3f;
          double bb_3f = g_1f3f * bb1f + g_2f3f * bb2f + g_3f3f * bb3f;

          // Calculate fluid-frame angle between photon direction and magnetic field
          double pf_dot_pf = p_1f * p1f + p_2f * p2f + p_3f * p3f;
          double bbf_dot_bbf = bb_1f * bb1f + bb_2f * bb2f + bb_3f * bb3f;
          double pf_dot_bbf = p_1f * bb1f + p_2f * bb2f + p_3f * bb3f;
          double cos2_theta = std::min(pf_dot_bbf * pf_dot_bbf / (pf_dot_pf * bbf_dot_bbf), 1.0);
          double sin_theta = std::sqrt(1.0 - cos2_theta);

          // Calculate fluid-frame quantities in CGS units
          double nu_fluid_cgs = (u0 * p_0 + u1 * p_1 + u2 * p_2 + u3 * p_3) / p_0 * im_freq;
          double n_cgs = rho * simulation_rho_cgs / (plasma_mu * physics::m_p);
          double n_e_cgs = n_cgs / (1.0 + 1.0 / plasma_ne_ni);
          double bb_cgs = std::sqrt(4.0 * math::pi * b_sq * e_unit);
          double nu_c_cgs = physics::e * bb_cgs / (2.0 * math::pi * physics::m_e * physics::c);
          double beta_inv = b_sq / (2.0 * pgas);
          double tt_rat = (plasma_rat_high + plasma_rat_low * beta_inv * beta_inv)
              / (1.0 + beta_inv * beta_inv);
          double kb_tt_tot_cgs = plasma_mu * physics::m_p * physics::c * physics::c * pgas / rho;
          double kb_tt_e_cgs = (plasma_ne_ni + 1.0) / (plasma_ne_ni + tt_rat) * kb_tt_tot_cgs;
          double theta_e = kb_tt_e_cgs / (physics::m_e * physics::c * physics::c);
          double nu_s_cgs = 2.0 / 9.0 * nu_c_cgs * theta_e * theta_e * sin_theta;

          // Calculate emission coefficient in CGS units (S 24), skipping contribution if necessary
          if (nu_s_cgs == 0.0)
            continue;
          double k2 = std::cyl_bessel_k(2.0, 1.0 / theta_e);
          if (k2 == 0.0)
            continue;
          double xx = nu_fluid_cgs / nu_s_cgs;
          double xx_factor = std::sqrt(xx) + std::sqrt(std::cbrt(32.0 * math::sqrt2 * xx));
          double j_nu_fluid_cgs = math::sqrt2 * math::pi * physics::e * physics::e * n_e_cgs
              * nu_s_cgs / (3.0 * k2 * physics::c) * xx_factor * xx_factor
              * std::exp(-std::cbrt(xx));
          if (j_nu_fluid_cgs == 0.0)
            continue;
          double nu_nu_fluid = nu_cgs / nu_fluid_cgs;
          double j_nu_cgs = nu_nu_fluid * nu_nu_fluid * j_nu_fluid_cgs;

          // Calculate absorption coefficient in CGS units (S 25)
          double b_nu_cgs = 2.0 * physics::h * nu_fluid_cgs * nu_fluid_cgs * nu_fluid_cgs
              / (physics::c * physics::c) / std::expm1(physics::h * nu_fluid_cgs / kb_tt_e_cgs);
          double k_nu_fluid_cgs = j_nu_fluid_cgs / b_nu_cgs;
          double k_nu_cgs = k_nu_fluid_cgs / nu_nu_fluid;

          // Calculate change in invariant intensity
          double i_nu_cgs = nu_cgs * nu_cgs * nu_cgs * image(m,l);
          double delta_s_cgs = delta_lambda * x_unit;
          double delta_tau_nu = k_nu_cgs * delta_s_cgs;
          double ss_nu_cgs = j_nu_cgs / k_nu_cgs;
          if (delta_tau_nu <= delta_tau_max)
            i_nu_cgs = std::exp(-delta_tau_nu) * (i_nu_cgs + ss_nu_cgs * std::expm1(delta_tau_nu));
          else
            i_nu_cgs = ss_nu_cgs;
          image(m,l) = i_nu_cgs / (nu_cgs * nu_cgs * nu_cgs);
        }
      }

    // Transform I_nu/nu^3 to brightness temperature
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        image(m,l) *= im_freq * physics::c * physics::c / (2.0 * physics::k_b);
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on formula
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_dir, and sample_len have been set.
//   Allocates and initializes image.
//   Assumes x^0 is ignorable.
//   References code comparison paper 2020 ApJ 897 148 (C).
void RayTracer::IntegrateFormulaRadiation()
{
  // Allocate image array
  image.Allocate(im_res, im_res);
  image.Zero();

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcon(4, 4);

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Go through samples
        int num_steps = sample_num(m,l);
        for (int n = 0; n < num_steps; n++)
        {
          // Check that this sample contributes
          double delta_lambda = sample_len(m,l,n);
          if (delta_lambda == 0.0)
            continue;

          // Extract geodesic position and momentum
          double x = sample_pos(m,l,n,1);
          double y = sample_pos(m,l,n,2);
          double z = sample_pos(m,l,n,3);
          double p_0 = sample_dir(m,l,n,0);
          double p_1 = sample_dir(m,l,n,1);
          double p_2 = sample_dir(m,l,n,2);
          double p_3 = sample_dir(m,l,n,3);

          // Calculate coordinates
          double r = RadialGeodesicCoordinate(x, y, z);
          double rr = std::sqrt(r * r - z * z);
          double cth = z / r;
          double sth = std::sqrt(1.0 - cth * cth);
          double ph = std::atan2(y, x) + std::atan(bh_a / r);
          double sph = std::sin(ph);
          double cph = std::cos(ph);

          // Calculate metric
          ContravariantGeodesicMetric(x, y, z, gcon);
          double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
          double sigma = r * r + bh_a * bh_a * cth * cth;
          double gtt_bl = -(1.0 + 2.0 * bh_m * r * (r * r + bh_a * bh_a) / (delta * sigma));
          double gtph_bl = -2.0 * bh_m * bh_a * r / (delta * sigma);
          double grr_bl = delta / sigma;
          double gthth_bl = 1.0 / sigma;
          double gphph_bl = (sigma - 2.0 * bh_m * r) / (delta * sigma * sth * sth);

          // Calculate frequency in CGS units
          double p0 = gcon(0,0) * p_0 + gcon(0,1) * p_1 + gcon(0,2) * p_2 + gcon(0,3) * p_3;
          double nu_cgs = -p0 / p_0 * im_freq;

          // Calculate angular momentum (C 6)
          double ll = formula_l0 / (1.0 + rr) * std::pow(rr, 1.0 + formula_q);

          // Calculate 4-velocity (C 7-8)
          double u_norm = 1.0 / std::sqrt(-gtt_bl + 2.0 * gtph_bl * ll - gphph_bl * ll * ll);
          double u_t_bl = -u_norm;
          double u_r_bl = 0.0;
          double u_th_bl = 0.0;
          double u_ph_bl = u_norm * ll;
          double ut_bl = gtt_bl * u_t_bl + gtph_bl * u_ph_bl;
          double ur_bl = grr_bl * u_r_bl;
          double uth_bl = gthth_bl * u_th_bl;
          double uph_bl = gtph_bl * u_t_bl + gphph_bl * u_ph_bl;
          double ut = ut_bl + 2.0 * bh_m * r / delta * ur_bl;
          double ur = ur_bl;
          double uth = uth_bl;
          double uph = uph_bl + bh_a / delta * ur_bl;
          double u0 = ut;
          double u1 = sth * cph * ur + cth * (r * cph + bh_a * sph) * uth
              + sth * (-r * sph + bh_a * cph) * uph;
          double u2 = sth * sph * ur + cth * (r * sph - bh_a * cph) * uth
              + sth * (r * cph + bh_a * sph) * uph;
          double u3 = cth * ur - r * sth * uth;

          // Calculate frequencies
          double nu_fluid_cgs = (u0 * p_0 + u1 * p_1 + u2 * p_2 + u3 * p_3) / p_0 * im_freq;
          double nu_nu_fluid = nu_cgs / nu_fluid_cgs;

          // Calculate fluid-frame number density (C 5)
          double n_n0_fluid = std::exp(-0.5
              * (r * r / (formula_r0 * formula_r0) + formula_h * formula_h * cth * cth));

          // Calculate emission coefficient in CGS units (C 9-10)
          double j_nu_fluid_cgs =
              formula_cn0 * n_n0_fluid * std::pow(nu_fluid_cgs / formula_nup, -formula_alpha);
          double j_nu_cgs = nu_nu_fluid * nu_nu_fluid * j_nu_fluid_cgs;

          // Calculate absorption coefficient in CGS units (C 11-12)
          double k_nu_fluid_cgs = formula_a * formula_cn0 * n_n0_fluid
              * std::pow(nu_fluid_cgs / formula_nup, -formula_beta - formula_alpha);
          double k_nu_cgs = k_nu_fluid_cgs / nu_nu_fluid;

          // Calculate change in invariant intensity
          double i_nu_cgs = nu_cgs * nu_cgs * nu_cgs * image(m,l);
          double delta_s_cgs = delta_lambda * formula_mass;
          if (k_nu_cgs > 0.0)
          {
            double delta_tau_nu = k_nu_cgs * delta_s_cgs;
            double ss_nu_cgs = j_nu_cgs / k_nu_cgs;
            if (delta_tau_nu <= delta_tau_max)
              i_nu_cgs = std::exp(-delta_tau_nu) * (i_nu_cgs + ss_nu_cgs * std::expm1(delta_tau_nu));
            else
              i_nu_cgs = ss_nu_cgs;
          } else
            i_nu_cgs += j_nu_cgs * delta_s_cgs;
          image(m,l) = i_nu_cgs / (nu_cgs * nu_cgs * nu_cgs);
        }
      }

    // Transform I_nu/nu^3 to brightness temperature
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        image(m,l) *= im_freq * physics::c * physics::c / (2.0 * physics::k_b);
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating radial coordinate given location in coordinates used for geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   returned value: radial coordinate
// Notes:
//   Assumes Cartesian Kerr-Schild coordinates.
double RayTracer::RadialGeodesicCoordinate(double x, double y, double z)
{
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
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
void RayTracer::CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov)
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
  double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

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
void RayTracer::ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon)
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
  double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

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
void RayTracer::ContravariantGeodesicMetricDerivative(double x, double y, double z,
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
  double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

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

//--------------------------------------------------------------------------------------------------

// Function for calculating covariant metric components in simulation coordinates
// Inputs:
//   x1, x2, x3: coordinates
// Output:
//   gcov: components set
// Notes:
//   Assumes gcov is allocated to be 4*4.
void RayTracer::CovariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcov)
{
  // Account for simulation metric
  switch (simulation_coord)
  {
    // Spherical Kerr-Schild
    case sph_ks:
    {
      // Calculate useful quantities
      double r = x1;
      double th = x2;
      double sth = std::sin(th);
      double cth = std::cos(th);
      double sigma = r * r + bh_a * bh_a * cth * cth;

      // Calculate metric components
      gcov(0,0) = -(1.0 - 2.0 * bh_m * r / sigma);
      gcov(0,1) = 2.0 * bh_m * r / sigma;
      gcov(0,2) = 0.0;
      gcov(0,3) = -2.0 * bh_m * bh_a * r * sth * sth / sigma;
      gcov(1,0) = 2.0 * bh_m * r / sigma;
      gcov(1,1) = 1.0 + 2.0 * bh_m * r / sigma;
      gcov(1,2) = 0.0;
      gcov(1,3) = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * sth * sth;
      gcov(2,0) = 0.0;
      gcov(2,1) = 0.0;
      gcov(2,2) = sigma;
      gcov(2,3) = 0.0;
      gcov(3,0) = -2.0 * bh_m * bh_a * r * sth * sth / sigma;
      gcov(3,1) = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * sth * sth;
      gcov(3,2) = 0.0;
      gcov(3,3) =
          (r * r + bh_a * bh_a + 2.0 * bh_m * bh_a * bh_a * r * sth * sth / sigma) * sth * sth;
      break;
    }

    // Cartesian Kerr-Schild
    case cart_ks:
    {
      // Calculate useful quantities
      double x = x1;
      double y = x2;
      double z = x3;
      double a2 = bh_a * bh_a;
      double rr2 = x * x + y * y + z * z;
      double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
      double r = std::sqrt(r2);
      double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

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
      break;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating contravariant metric components in simulation coordinates
// Inputs:
//   x1, x2, x3: coordinates
// Output:
//   gcon: components set
// Notes:
//   Assumes gcon is allocated to be 4*4.
void RayTracer::ContravariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcon)
{
  // Account for simulation metric
  switch (simulation_coord)
  {
    // Spherical Kerr-Schild
    case sph_ks:
    {
      // Calculate useful quantities
      double r = x1;
      double th = x2;
      double sth = std::sin(th);
      double cth = std::cos(th);
      double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
      double sigma = r * r + bh_a * bh_a * cth * cth;

      // Calculate metric components
      gcon(0,0) = -(1.0 + 2.0 * bh_m * r / sigma);
      gcon(0,1) = 2.0 * bh_m * r / sigma;
      gcon(0,2) = 0.0;
      gcon(0,3) = 0.0;
      gcon(1,0) = 2.0 * bh_m * r / sigma;
      gcon(1,1) = delta / sigma;
      gcon(1,2) = 0.0;
      gcon(1,3) = bh_a / sigma;
      gcon(2,0) = 0.0;
      gcon(2,1) = 0.0;
      gcon(2,2) = 1.0 / sigma;
      gcon(2,3) = 0.0;
      gcon(3,0) = 0.0;
      gcon(3,1) = bh_a / sigma;
      gcon(3,2) = 0.0;
      gcon(3,3) = 1.0 / (sigma * sth * sth);
      break;
    }

    // Cartesian Kerr-Schild
    case cart_ks:
    {
      // Calculate useful quantities
      double x = x1;
      double y = x2;
      double z = x3;
      double a2 = bh_a * bh_a;
      double rr2 = x * x + y * y + z * z;
      double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
      double r = std::sqrt(r2);
      double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

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
      break;
    }
  }
  return;
}
