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
  double im_r;
  double im_th;
  double im_ph;
  double im_rot;
  double im_width;
  int im_res;
  double im_step;
  int im_max_steps;

  // Grid data
  double r_min, r_max, th_min, th_max, ph_min, ph_max;
  array<float> rf, thf, phf;
  array<float> rho;

  // Sample and image data
  double r_hor;
  const int ind_rho = 0;
  const int num_states = 1;
  array<double> im_pos, im_dir;
  array<double> sample_pos, sample_dir, sample_len;

  // Functions
  void make_image();
  void initialize_camera();
  void initialize_geodesics();
  void integrate_geodesics();
  void gcov_func(double r, double th, array<double> &gcov);
  void gcon_func(double r, double th, array<double> &gcon);
  void dgcon_func(double r, double th, array<double> &dgcon);
};

#endif
