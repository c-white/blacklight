// Ray Trace ray tracer header

#ifndef RAY_TRACER_H_
#define RAY_TRACER_H_

// Ray Trace headers
#include "array.hpp"        // array
#include "read_athena.hpp"  // athena_reader
#include "read_input.hpp"   // input_reader

//--------------------------------------------------------------------------------------------------

// Ray tracer
struct ray_tracer
{
  // Constructors and destructor
  ray_tracer(const input_reader &inputs, const athena_reader &raw_data);
  ray_tracer(const ray_tracer &source) = delete;
  ray_tracer &operator=(const ray_tracer &source) = delete;
  ~ray_tracer();

  // Input data - coordinates
  double m;
  double a;

  // Input data - image
  double im_th;
  double im_ph;
  double im_width;
  int im_res;
  int num_samples;

  // Grid data
  array<float> rf, thf, phf;
  array<float> rho;

  // Sample and image data
  array<float> samples;
  array<float> image;

  // Functions
  void sample_rays();
  void make_image();
};

#endif
