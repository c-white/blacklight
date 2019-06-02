// Ray Trace ray tracer

// C++ headers
#include <cmath>  // acos, atan2, cos, hypot, sin

// Ray Trace headers
#include "ray_tracer.hpp"
#include "array.hpp"        // array
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
  im_th = inputs.im_th;
  im_ph = inputs.im_ph;
  im_rot = inputs.im_rot;
  im_width = inputs.im_width;
  im_res = inputs.im_res;
  num_samples = inputs.num_samples;

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
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for setting up camera pixels and initial ray directions
// Inputs: (none)
// Output: (none)
// Notes:
//   Allocates and initializes im_pos and im_dir.
//   Neglects spacetime curvature at camera location.
//   Symbols:
//     n: unit outward normal
//     u: unit right vector
//     v: unit up vector
void ray_tracer::initialize_camera()
{
  // Calculate camera position
  double im_r = r_max;
  double im_x = r_max * std::sin(im_th) * std::cos(im_ph);
  double im_y = r_max * std::sin(im_th) * std::sin(im_ph);
  double im_z = r_max * std::cos(im_th);

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
  im_pos.allocate(3, im_res, im_res);
  im_dir.allocate(3, im_res, im_res);

  // Initialize arrays
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
    {
      // Calculate position in camera plane
      double u = (l - im_res/2.0 + 0.5) * bh_m * im_width / im_res;
      double v = (m - im_res/2.0 + 0.5) * bh_m * im_width / im_res;

      // Calculate position in space
      double x = im_x + u * im_ux + v * im_vx;
      double y = im_y + u * im_uy + v * im_vy;
      double z = im_z + u * im_uz + v * im_vz;
      double r = std::hypot(x, y, z);
      double th = std::acos(z / r);
      double ph = std::atan2(y, x);
      im_pos(0,m,l) = r;
      im_pos(1,m,l) = th;
      im_pos(2,m,l) = ph;

      // Calculate direction
      double nr = (x * im_nx + y * im_ny + z * im_nz) / r;
      double nth = (std::cos(th) * std::cos(ph) * im_nx + std::cos(th) * std::sin(ph) * im_ny
          - std::sin(th) * im_nz) / r;
      double nph = (-std::sin(ph) / std::sin(th) * im_nx + std::cos(ph) / std::sin(th) * im_ny) / r;
      im_dir(0,m,l) = nr;
      im_dir(1,m,l) = nth;
      im_dir(2,m,l) = nph;
    }
  return;
}
