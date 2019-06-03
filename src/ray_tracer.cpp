// Ray Trace ray tracer

// C++ headers
#include <cmath>  // acos, atan2, copysign, cos, fmax, hypot, sin, sqrt

// Ray Trace headers
#include "ray_tracer.hpp"
#include "array.hpp"        // array
#include "exceptions.hpp"   // ray_trace_exception
#include "read_athena.hpp"  // athena_reader
#include "read_input.hpp"   // input_reader

//--------------------------------------------------------------------------------------------------

// Ray tracer constructor
// Inputs:
//   inputs: object containing input parameters read from input file
//   raw_data: object containing raw data read from data file
ray_tracer::ray_tracer(const input_reader &inputs, const athena_reader &raw_data)
{
  // Copy coordinate input data
  bh_m = inputs.bh_m;
  bh_a = inputs.bh_a;

  // Copy image input data
  im_r = inputs.im_r;
  im_th = inputs.im_th;
  im_ph = inputs.im_ph;
  im_rot = inputs.im_rot;
  im_width = inputs.im_width;
  im_res = inputs.im_res;
  im_step = inputs.im_step;
  im_max_steps = inputs.im_max_steps;

  // Copy raw data scalars
  r_min = raw_data.r_min;
  r_max = raw_data.r_max;
  th_min = raw_data.th_min;
  th_max = raw_data.th_max;
  ph_min = raw_data.ph_min;
  ph_max = raw_data.ph_max;

  // Make shallow copies of raw data arrays
  rf = raw_data.rf;
  thf = raw_data.thf;
  phf = raw_data.phf;
  rho = raw_data.prim;
  rho.slice(5, raw_data.ind_rho);

  // Calculate horizon radius
  r_hor = bh_m + std::sqrt(bh_m * bh_m - bh_a * bh_a);
}

//--------------------------------------------------------------------------------------------------

// Ray tracer destructor
ray_tracer::~ray_tracer() {}

//--------------------------------------------------------------------------------------------------

// Top-level function for processing raw data into image
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes all data arrays have been set.
void ray_tracer::make_image()
{
  initialize_camera();
  initialize_geodesics();
  integrate_geodesics();
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
void ray_tracer::initialize_camera()
{
  // Calculate camera position
  double im_x = im_r * std::sin(im_th) * std::cos(im_ph);
  double im_y = im_r * std::sin(im_th) * std::sin(im_ph);
  double im_z = im_r * std::cos(im_th);

  // Calculate camera direction
  double im_nx = im_x / im_r;
  double im_ny = im_y / im_r;
  double im_nz = im_z / im_r;

  // Calculate camera orientation
  double im_uthh = std::sin(im_rot);
  double im_uphh = std::cos(im_rot);
  double im_ux =
      std::cos(im_th) * std::cos(im_ph) * im_uthh - std::sin(im_th) * std::sin(im_ph) * im_uphh;
  double im_uy =
      std::cos(im_th) * std::sin(im_ph) * im_uthh + std::sin(im_th) * std::cos(im_ph) * im_uphh;
  double im_uz = -std::sin(im_th) * im_uthh;
  double im_vthh = -std::cos(im_rot);
  double im_vphh = std::sin(im_rot);
  double im_vx =
      std::cos(im_th) * std::cos(im_ph) * im_vthh - std::sin(im_th) * std::sin(im_ph) * im_vphh;
  double im_vy =
      std::cos(im_th) * std::sin(im_ph) * im_vthh + std::sin(im_th) * std::cos(im_ph) * im_vphh;
  double im_vz = -std::sin(im_th) * im_vthh;

  // Allocate arrays
  im_pos.allocate(im_res, im_res, 4);
  im_dir.allocate(im_res, im_res, 4);

  // Initialize arrays
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
    {
      // Calculate position in camera plane
      double u = (l - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
      double v = (m - im_res/2.0 + 0.5) * bh_m * im_width / im_res;

      // Calculate position in space
      double t = 0.0;
      double x = im_x + u * im_ux + v * im_vx;
      double y = im_y + u * im_uy + v * im_vy;
      double z = im_z + u * im_uz + v * im_vz;
      double r = std::hypot(x, y, z);
      double th = std::acos(z / r);
      double ph = std::atan2(y, x);
      im_pos(m,l,0) = t;
      im_pos(m,l,1) = r;
      im_pos(m,l,2) = th;
      im_pos(m,l,3) = ph;

      // Calculate direction
      double nr = (x * im_nx + y * im_ny + z * im_nz) / r;
      double nth = (std::cos(th) * std::cos(ph) * im_nx + std::cos(th) * std::sin(ph) * im_ny
          - std::sin(th) * im_nz) / r;
      double nph = (-std::sin(ph) / std::sin(th) * im_nx + std::cos(ph) / std::sin(th) * im_ny) / r;
      im_dir(m,l,1) = nr;
      im_dir(m,l,2) = nth;
      im_dir(m,l,3) = nph;
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
void ray_tracer::initialize_geodesics()
{
  // Allocate scratch array
  array<double> gcov(4, 4);

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
      gcov_func(x[1], x[2], gcov);
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
//   Allocates and initializes sample_pos, sample_dir, and sample_len.
//   Assumes x^0 and x^3 are ignorable coordinates.
//   Integrates via the midpoint method (2nd-order RK).
void ray_tracer::integrate_geodesics()
{
  // Allocate arrays
  sample_pos.allocate(im_res, im_res, im_max_steps, 4);
  sample_dir.allocate(im_res, im_res, im_max_steps, 4);
  sample_len.allocate(im_res, im_res, im_max_steps);
  sample_len.zero();

  // Allocate scratch arrays
  array<double> gcon(4, 4);
  array<double> dgcon(2, 4, 4);

  // Go through image pixels
  bool max_steps_exceeded = false;
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
        double step = -im_step * x[1];

        // Calculate position at half step, checking that step is worth taking
        gcon_func(x[1], x[2], gcon);
        double dx1[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            dx1[mu] += gcon(mu,nu) * p[nu];
        if ((x[1] > im_r and dx1[1] < 0.0) or x[1] < r_hor)
          break;
        for (int mu = 0; mu < 4; mu++)
          sample_pos(m,l,n,mu) = x[mu] + step/2.0 * dx1[mu];

        // Calculate momentum at half step
        dgcon_func(x[1], x[2], dgcon);
        double dp1[4] = {};
        for (int a = 1; a <= 2; a++)
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              dp1[a] -= 0.5 * dgcon(a-1,mu,nu) * p[mu] * p[nu];
        for (int mu = 0; mu < 4; mu++)
          sample_dir(m,l,n,mu) = p[mu] + step/2.0 * dp1[mu];

        // Calculate position at full step
        gcon_func(sample_pos(m,l,n,1), sample_pos(m,l,n,2), gcon);
        double dx2[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            dx2[mu] += gcon(mu,nu) * sample_dir(m,l,n,nu);
        for (int mu = 0; mu < 4; mu++)
          x[mu] += step * dx2[mu];

        // Calculate momentum at full step
        dgcon_func(sample_pos(m,l,n,1), sample_pos(m,l,n,2), dgcon);
        double dp2[4] = {};
        for (int a = 1; a <= 2; a++)
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              dp2[a] -= 0.5 * dgcon(a-1,mu,nu) * sample_dir(m,l,n,mu) * sample_dir(m,l,n,nu);
        for (int mu = 0; mu < 4; mu++)
          p[mu] += step * dp2[mu];

        // Store length of step
        sample_len(m,l,n) = -step;

        // Check for too many steps taken
        if (n == im_max_steps - 1 and not ((x[1] > im_r and dx1[1] < 0.0) or x[1] < r_hor))
          max_steps_exceeded = true;
      }
    }
  if (max_steps_exceeded)
    throw ray_trace_exception("Error: Geodesic requires too many steps\n");
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating covariant metric components
// Inputs:
//   r: radial coordinate
//   th: polar coordinate
// Output:
//   gcov: components set
// Notes:
//   Assumes gcov is allocated to be 4*4.
//   Assumes spherical Kerr-Schild coordinates.
void ray_tracer::gcov_func(double r, double th, array<double> &gcov)
{
  double sth = std::sin(th);
  double s2th = sth * sth;
  double c2th = 1.0 - s2th;
  double sigma = r * r + bh_a * bh_a * c2th;
  gcov(0,0) = -(1.0 - 2.0 * bh_m * r / sigma);
  gcov(0,1) = gcov(1,0) = 2.0 * bh_m * r / sigma;
  gcov(0,2) = gcov(2,0) = 0.0;
  gcov(0,3) = gcov(3,0) = -2.0 * bh_m * bh_a * r * s2th / sigma;
  gcov(1,1) = 1.0 + 2.0 * bh_m * r / sigma;
  gcov(1,2) = gcov(2,1) = 0.0;
  gcov(1,3) = gcov(3,1) = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * s2th;
  gcov(2,2) = sigma;
  gcov(2,3) = gcov(3,2) = 0.0;
  gcov(3,3) = (r * r + bh_a * bh_a + 2.0 * bh_m * bh_a * bh_a * r * s2th / sigma) * s2th;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating contravariant metric components
// Inputs:
//   r: radial coordinate
//   th: polar coordinate
// Output:
//   gcon: components set
// Notes:
//   Assumes gcon is allocated to be 4*4.
//   Assumes spherical Kerr-Schild coordinates.
void ray_tracer::gcon_func(double r, double th, array<double> &gcon)
{
  double sth = std::sin(th);
  double s2th = sth * sth;
  double c2th = 1.0 - s2th;
  double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
  double sigma = r * r + bh_a * bh_a * c2th;
  gcon(0,0) = -(1.0 + 2.0 * bh_m * r / sigma);
  gcon(0,1) = gcon(1,0) = 2.0 * bh_m * r / sigma;
  gcon(0,2) = gcon(2,0) = 0.0;
  gcon(0,3) = gcon(3,0) = 0.0;
  gcon(1,1) = delta / sigma;
  gcon(1,2) = gcon(2,1) = 0.0;
  gcon(1,3) = gcon(3,1) = bh_a / sigma;
  gcon(2,2) = 1.0 / sigma;
  gcon(2,3) = gcon(3,2) = 0.0;
  gcon(3,3) = 1.0 / (sigma * s2th);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating derivatives of contravariant metric components
// Inputs:
//   r: radial coordinate
//   th: polar coordinate
// Output:
//   dgcon: components set
// Notes:
//   Assumes dgcon is allocated to be 2*4*4.
//   Assumes spherical Kerr-Schild coordinates.
void ray_tracer::dgcon_func(double r, double th, array<double> &dgcon)
{
  double sth = std::sin(th);
  double cth = std::cos(th);
  double s2th = sth * sth;
  double c2th = cth * cth;
  double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
  double sigma = r * r + bh_a * bh_a * c2th;
  double sigma2 = sigma * sigma;
  dgcon(0,0,0) = -2.0 * bh_m * (sigma - 2.0 * r * r) / sigma2;
  dgcon(0,0,1) = dgcon(0,1,0) = 2.0 * bh_m * (sigma - 2.0 * r * r) / sigma2;
  dgcon(0,0,2) = dgcon(0,2,0) = 0.0;
  dgcon(0,0,3) = dgcon(0,3,0) = 0.0;
  dgcon(0,1,1) = 2.0 * ((r - bh_m) * sigma - r * delta) / sigma2;
  dgcon(0,1,2) = dgcon(0,2,1) = 0.0;
  dgcon(0,1,3) = dgcon(0,3,1) = -2.0 * bh_a * r / sigma2;
  dgcon(0,2,2) = -2.0 * r / sigma2;
  dgcon(0,2,3) = dgcon(0,3,2) = 0.0;
  dgcon(0,3,3) = -2.0 * r / (sigma2 * s2th);
  dgcon(1,0,0) = -4.0 * bh_m * bh_a * bh_a * r * sth * cth / sigma2;
  dgcon(1,0,1) = dgcon(1,1,0) = 4.0 * bh_m * bh_a * bh_a * r * sth * cth / sigma2;
  dgcon(1,0,2) = dgcon(1,2,0) = 0.0;
  dgcon(1,0,3) = dgcon(1,3,0) = 0.0;
  dgcon(1,1,1) = 2.0 * bh_a * bh_a * sth * cth * delta / sigma2;
  dgcon(1,1,2) = dgcon(1,2,1) = 0.0;
  dgcon(1,1,3) = dgcon(1,3,1) = 2.0 * bh_a * bh_a * bh_a * sth * cth / sigma2;
  dgcon(1,2,2) = 2.0 * bh_a * bh_a * sth * cth / sigma2;
  dgcon(1,2,3) = dgcon(1,3,2) = 0.0;
  dgcon(1,3,3) = 2.0 * cth * (bh_a * bh_a * s2th - sigma) / (sigma2 * sth * s2th);
  return;
}
