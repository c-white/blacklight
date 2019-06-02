// Ray Trace ray tracer header

#ifndef RAY_TRACER_H_
#define RAY_TRACER_H_

// Ray Trace headers
#include "array.hpp"        // array
#include "read_athena.hpp"  // athena_reader

//--------------------------------------------------------------------------------------------------

// Ray tracer
struct ray_tracer
{
  // Constructors and destructor
  ray_tracer(const athena_reader &raw_data);
  ray_tracer(const ray_tracer &source) = delete;
  ray_tracer &operator=(const ray_tracer &source) = delete;
  ~ray_tracer();

  // Data
  array<float> rf, thf, phf;
  array<float> rho;

  // Functions
};

#endif
