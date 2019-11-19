// Ray Trace ray tracer

// C++ headers
#include <algorithm>  // max
#include <cmath>      // acos, atan, atan2, copysign, cos, fmax, fmod, hypot, sin, sqrt
#include <sstream>    // stringstream
#include <string>     // string

// Ray Trace headers
#include "ray_tracer.hpp"
#include "array.hpp"        // Array
#include "exceptions.hpp"   // RayTraceWarning
#include "ray_trace.hpp"    // math
#include "read_athena.hpp"  // AthenaReader
#include "read_input.hpp"   // InputReader

//--------------------------------------------------------------------------------------------------

// Ray tracer constructor
// Inputs:
//   input_reader: object containing input parameters read from input file
//   athena_reader: object containing raw data read from data file
RayTracer::RayTracer(const InputReader &input_reader, const AthenaReader &athena_reader)
{
  // Copy coordinate input data
  bh_m = input_reader.bh_m;
  bh_a = input_reader.bh_a;

  // Copy image input data
  im_r = input_reader.im_r;
  im_th = input_reader.im_th;
  im_ph = input_reader.im_ph;
  im_rot = input_reader.im_rot;
  im_width = input_reader.im_width;
  im_res = input_reader.im_res;
  im_step = input_reader.im_step;
  im_max_steps = input_reader.im_max_steps;

  // Copy raw data scalars
  r_min = athena_reader.r_min;
  r_max = athena_reader.r_max;
  th_min = athena_reader.th_min;
  th_max = athena_reader.th_max;
  ph_min = athena_reader.ph_min;
  ph_max = athena_reader.ph_max;

  // Make shallow copies of raw data arrays
  rf = athena_reader.rf;
  thf = athena_reader.thf;
  phf = athena_reader.phf;
  rho = athena_reader.prim;
  rho.Slice(5, athena_reader.ind_rho);

  // Calculate horizon radius
  r_hor = bh_m + std::sqrt(bh_m * bh_m - bh_a * bh_a);
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
  SampleAlongGeodesics();
  IntegrateRadiation();
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
  double im_x = im_sth * (im_r * im_cph - bh_a * im_sph);
  double im_y = im_sth * (im_r * im_sph + bh_a * im_cph);
  double im_z = im_r * im_cth;

  // Calculate camera direction
  double im_nx = im_sth * im_cph;
  double im_ny = im_sth * im_sph;
  double im_nz = im_cth;

  // Calculate camera orientation
  double im_uthh = im_srot;
  double im_uphh = im_crot;
  double im_ux = im_cth * im_cph * im_uthh - im_sth * im_sph * im_uphh;
  double im_uy = im_cth * im_sph * im_uthh + im_sth * im_cph * im_uphh;
  double im_uz = -im_sth * im_uthh;
  double im_vthh = -im_crot;
  double im_vphh = im_srot;
  double im_vx = im_cth * im_cph * im_vthh - im_sth * im_sph * im_vphh;
  double im_vy = im_cth * im_sph * im_vthh + im_sth * im_cph * im_vphh;
  double im_vz = -im_sth * im_vthh;

  // Allocate arrays
  im_pos.Allocate(im_res, im_res, 4);
  im_dir.Allocate(im_res, im_res, 4);

  // Initialize arrays
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
    {
      // Calculate position in camera plane
      double u = (l - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
      double v = (m - im_res/2.0 + 0.5) * bh_m * im_width / im_res;

      // Calculate position in space
      im_pos(m,l,0) = 0.0;
      im_pos(m,l,1) = im_x + u * im_ux + v * im_vx;
      im_pos(m,l,2) = im_y + u * im_uy + v * im_vy;
      im_pos(m,l,3) = im_z + u * im_uz + v * im_vz;

      // Set direction
      im_dir(m,l,1) = im_nx;
      im_dir(m,l,2) = im_ny;
      im_dir(m,l,3) = im_nz;
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
  // Allocate scratch array
  Array<double> gcov(4, 4);

  // Go through image pixels
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
      for (int mu = 1; mu < 4; mu++)
      {
        im_dir(m,l,mu) = 0.0;
        for (int nu = 1; nu < 4; nu++)
          im_dir(m,l,mu) += gcov(mu,nu) * p[nu];
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
//   Allocates and initializes sample_pos, sample_dir, sample_len, and geodesic_flags.
//   Assumes x^0 is ignorable.
//   Integrates via the midpoint method (2nd-order RK).
// TODO: calculate better step size
void RayTracer::IntegrateGeodesics()
{
  // Allocate arrays
  sample_pos.Allocate(im_res, im_res, im_max_steps, 4);
  sample_dir.Allocate(im_res, im_res, im_max_steps, 4);
  sample_len.Allocate(im_res, im_res, im_max_steps);
  sample_len.Zero();
  geodesic_flags.Allocate(im_res, im_res);
  geodesic_flags.Zero();

  // Allocate scratch arrays
  Array<double> gcon(4, 4);
  Array<double> dgcon(3, 4, 4);

  // Go through image pixels
  im_steps = 0;
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
      for (int n = 0; n < im_max_steps; n++)
      {
        // Calculate step size for going back to source
        double r = RadialGeodesicCoordinate(x[1], x[2], x[3]);
        double step = -im_step * r;

        // Calculate position at half step, checking that step is worth taking
        ContravariantGeodesicMetric(x[1], x[2], x[3], gcon);
        double dx1[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            dx1[mu] += gcon(mu,nu) * p[nu];
        for (int mu = 0; mu < 4; mu++)
          sample_pos(m,l,n,mu) = x[mu] + step/2.0 * dx1[mu];
        double delta_r = RadialGeodesicCoordinate(sample_pos(m,l,n,1), sample_pos(m,l,n,2),
            sample_pos(m,l,n,3)) - r;
        if ((r > im_r and delta_r > 0.0) or r < r_hor)
          break;

        // Calculate momentum at half step
        ContravariantGeodesicMetricDerivative(x[1], x[2], x[3], dgcon);
        double dp1[4] = {};
        for (int a = 1; a <= 3; a++)
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              dp1[a] -= 0.5 * dgcon(a-1,mu,nu) * p[mu] * p[nu];
        for (int mu = 0; mu < 4; mu++)
          sample_dir(m,l,n,mu) = p[mu] + step/2.0 * dp1[mu];

        // Calculate position at full step
        ContravariantGeodesicMetric(sample_pos(m,l,n,1), sample_pos(m,l,n,2), sample_pos(m,l,n,3),
            gcon);
        double dx2[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            dx2[mu] += gcon(mu,nu) * sample_dir(m,l,n,nu);
        for (int mu = 0; mu < 4; mu++)
          x[mu] += step * dx2[mu];

        // Calculate momentum at full step
        ContravariantGeodesicMetricDerivative(sample_pos(m,l,n,1), sample_pos(m,l,n,2),
            sample_pos(m,l,n,3), dgcon);
        double dp2[4] = {};
        for (int a = 1; a <= 3; a++)
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              dp2[a] -= 0.5 * dgcon(a-1,mu,nu) * sample_dir(m,l,n,mu) * sample_dir(m,l,n,nu);
        for (int mu = 0; mu < 4; mu++)
          p[mu] += step * dp2[mu];

        // Store length of step
        sample_len(m,l,n) = -step;

        // Check for too many steps taken
        if (n == im_max_steps - 1 and not ((x[1] > im_r and dx1[1] < 0.0) or x[1] < r_hor))
          geodesic_flags(m,l) = true;
        im_steps = std::max(im_steps, n + 1);
      }
    }

  // Note how many geodesics do not terminate properly
  int num_bad_geodesics = 0;
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
      if (geodesic_flags(m,l))
        num_bad_geodesics++;
  if (num_bad_geodesics > 0)
  {
    std::stringstream message;
    message << num_bad_geodesics << " out of " << im_res * im_res
        << " geodesics terminate unexpectedly.";
    RayTraceWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for transforming geodesics from integrating metric to simulation metric
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_steps, sample_pos, sample_dir, and sample_len have been set.
//   Transforms sample_pos, sample_dir, and sample_len from integrating metric to simulation metric.
//   Assumes integrating metric is Cartesian Kerr-Schild.
//   Assumes simulation metric is spherical Kerr-Schild.
//   Transformation of time components is trivial.
//   Transformation of sample_len is trivial.
void RayTracer::TransformGeodesics()
{
  // Go through samples
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
      for (int n = 0; n < im_steps; n++)
      {
        // Extract Cartesian position
        double x = sample_pos(m,l,n,1);
        double y = sample_pos(m,l,n,2);
        double z = sample_pos(m,l,n,3);

        // Extract Cartesian direction
        double ux = sample_dir(m,l,n,1);
        double uy = sample_dir(m,l,n,2);
        double uz = sample_dir(m,l,n,3);

        // Calculate spherical position
        double a2 = bh_a * bh_a;
        double rr2 = x * x + y * y + z * z;
        double r2 = 0.5 * (rr2 - a2 + std::sqrt((rr2 - a2) * (rr2 - a2) + 4.0 * a2 * z * z));
        double r = std::sqrt(r2);
        double th = std::acos(z / r);
        double ph = std::atan2(y, x) - std::atan(bh_a / r);
        sample_pos(m,l,n,1) = r;
        sample_pos(m,l,n,2) = th;
        sample_pos(m,l,n,3) = ph;

        // Calculate Jacobian of transformation
        double dr_dx = r * x / (2.0 * r2 - rr2 + a2);
        double dr_dy = r * y / (2.0 * r2 - rr2 + a2);
        double dr_dz = (r * z + a2 * z / r) / (2.0 * r2 - rr2 + a2);
        double dth_dx = z * dr_dx / (r * std::sqrt(r2 - z * z));
        double dth_dy = z * dr_dy / (r * std::sqrt(r2 - z * z));
        double dth_dz = (z * dr_dz - r) / (r * std::sqrt(r2 - z * z));
        double dph_dx = bh_a * dr_dx - y / ((r2 + a2) * (1.0 - z * z / r2));
        double dph_dy = bh_a * dr_dy + x / ((r2 + a2) * (1.0 - z * z / r2));
        double dph_dz = bh_a * dr_dz;

        // Calculate spherical direction
        sample_dir(m,l,n,1) = dr_dx * ux + dr_dy * uy + dr_dz * uz;
        sample_dir(m,l,n,2) = dth_dx * ux + dth_dy * uy + dth_dz * uz;
        sample_dir(m,l,n,3) = dph_dx * ux + dph_dy * uy + dph_dz * uz;
      }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for resampling cell data onto rays.
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_steps, sample_len, and geodesic_flags have been set.
//   Assumes sample_pos has been transformed.
//   Allocates and initializes sample_rho.
void RayTracer::SampleAlongGeodesics()
{
  // Allocate resampling arrays
  sample_rho.Allocate(im_res, im_res, im_steps);

  // Prepare bookkeeping
  int n_b = rf.n2;
  int n_i = rf.n1 - 1;
  int n_j = thf.n1 - 1;
  int n_k = phf.n1 - 1;
  int b = 0;
  int i = 0;
  int j = 0;
  int k = 0;
  double r_min_block = rf(b,0);
  double r_max_block = rf(b,n_i);
  double th_min_block = thf(b,0);
  double th_max_block = thf(b,n_j);
  double ph_min_block = phf(b,0);
  double ph_max_block = phf(b,n_k);

  // Resample cell data onto geodesics
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
    {
      // Set fallback values if geodesic poorly terminated
      if (geodesic_flags(m,l))
      {
        for (int n = 0; n < im_steps; n++)
          sample_rho(m,l,n) = rho_fallback;
        continue;
      }

      // Go along geodesic
      for (int n = 0; n < im_steps; n++)
      {
        // End if geodesic terminated
        if (sample_len(m,l,n) == 0.0)
          continue;

        // Extract coordinates
        double r = sample_pos(m,l,n,1);
        double th = sample_pos(m,l,n,2);
        double ph = sample_pos(m,l,n,3);

        // Determine block
        if (r < r_min_block or r > r_max_block or th < th_min_block or th > th_max_block
            or ph < ph_min_block or ph > ph_max_block)
        {
          // Check if block contains position
          for (b = 0; b < n_b; b++)
          {
            double r_min_temp = rf(b,0);
            double r_max_temp = rf(b,n_i);
            if (r < r_min_temp or r > r_max_temp)
              continue;
            double th_min_temp = thf(b,0);
            double th_max_temp = thf(b,n_j);
            if (th < th_min_temp or th > th_max_temp)
              continue;
            double ph_min_temp = phf(b,0);
            double ph_max_temp = phf(b,n_k);
            if (ph < ph_min_temp or ph > ph_max_temp)
              continue;
            i = 0;
            j = 0;
            k = 0;
            r_min_block = r_min_temp;
            r_max_block = r_max_temp;
            th_min_block = th_min_temp;
            th_max_block = th_max_temp;
            ph_min_block = ph_min_temp;
            ph_max_block = ph_max_temp;
          }

          // Set fallback values if off grid
          if (b == n_b)
          {
            sample_rho(m,l,n) = rho_fallback;
            continue;
          }
        }

        // Determine cell
        for (i = 0; i < n_i; i++)
          if (static_cast<double>(rf(b,i+1)) >= r)
            break;
        for (j = 0; j < n_j; j++)
          if (static_cast<double>(thf(b,j+1)) >= th)
            break;
        for (k = 0; k < n_k; k++)
          if (static_cast<double>(phf(b,k+1)) >= ph)
            break;

        // Resample values
        sample_rho(m,l,n) = rho(b,k,j,i);
      }
    }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_steps and sample_rho have been set.
//   Assumes sample_dir and sample_len have been transformed.
//   Allocates and initializes image.
//   TODO: use physically meaningful formula
void RayTracer::IntegrateRadiation()
{
  // Allocate array
  image.Allocate(im_res, im_res);
  image.Zero();

  // Integrate radiative transfer equation
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
      for (int n = 0; n < im_steps; n++)
        image(m,l) += sample_rho(m,l,n) * static_cast<float>(sample_len(m,l,n));
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
//   Assumes Cartesian Kerr-Schild coordinates.
void RayTracer::CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov)
{
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
//   Assumes Cartesian Kerr-Schild coordinates.
void RayTracer::ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon)
{
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
//   Assumes Cartesian Kerr-Schild coordinates.
void RayTracer::ContravariantGeodesicMetricDerivative(double x, double y, double z,
    Array<double> &dgcon)
{
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
