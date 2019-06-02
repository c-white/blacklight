// Ray Trace ray tracer

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
  m = inputs.m;
  a = inputs.a;

  // Copy image input data
  im_th = inputs.im_th;
  im_ph = inputs.im_ph;
  im_width = inputs.im_width;
  im_res = inputs.im_res;
  num_samples = inputs.num_samples;

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

// Top-level function for integrating geodesics and extracting grid information
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes all data arrays have been set.
//   Allocates and initializes samples.
void ray_tracer::sample_rays()
{
  return;
}

//--------------------------------------------------------------------------------------------------

// Top-level function for creating image from sampled geodesics
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes all samples has been set.
//   Allocates and initializes image.
void ray_tracer::make_image()
{
  return;
}
