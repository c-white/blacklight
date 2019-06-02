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
  double bh_m;
  double bh_a;

  // Input data - image
  double im_th;
  double im_ph;
  double im_rot;
  double im_width;
  int im_res;
  int num_samples;

  // Grid data
  double r_min, r_max, th_min, th_max, ph_min, ph_max;
  array<float> rf, thf, phf;
  array<float> rho;

  // Sample and image data
  const int ind_rho = 0;
  const int num_states = 1;
  array<double> im_pos, im_dir;

  // Functions
  void make_image();
  void initialize_camera();
};

#endif
