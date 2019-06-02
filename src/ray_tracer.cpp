// Ray Trace ray tracer

// Ray Trace headers
#include "ray_tracer.hpp"
#include "array.hpp"        // array
#include "read_athena.hpp"  // athena_reader

//--------------------------------------------------------------------------------------------------

// Ray tracer constructor
// Inputs:
//   raw_data: object containing raw data read from data file
ray_tracer::ray_tracer(const athena_reader &raw_data)
{
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
