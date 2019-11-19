// Ray Trace ray tracer header

#ifndef RAY_TRACER_H_
#define RAY_TRACER_H_

// Ray Trace headers
#include "array.hpp"        // Array
#include "read_athena.hpp"  // AthenaReader
#include "read_input.hpp"   // InputReader

//--------------------------------------------------------------------------------------------------

// Ray tracer
struct RayTracer
{
  // Constructors and destructor
  RayTracer(const InputReader &input_reader, const AthenaReader &athena_reader);
  RayTracer(const RayTracer &source) = delete;
  RayTracer &operator=(const RayTracer &source) = delete;
  ~RayTracer();

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

  // Input data - rays
  double ray_step;
  int ray_max_steps;
  bool flat;

  // Grid data
  double r_min, r_max, th_min, th_max, ph_min, ph_max;
  Array<float> rf, thf, phf;
  Array<float> rho;
  Array<float> pgas;

  // Fallback data
  const float rho_fallback = 1.0e-6f;
  const float pgas_fallback = 1.0e-8f;

  // Sample and image data
  double r_hor;
  int im_steps;
  Array<double> im_pos, im_dir;
  Array<double> sample_pos, sample_dir, sample_len;
  Array<bool> geodesic_flags;
  Array<float> sample_rho;
  Array<float> sample_pgas;
  Array<float> image;

  // Functions
  void MakeImage();
  void InitializeCamera();
  void InitializeGeodesics();
  void IntegrateGeodesics();
  void TransformGeodesics();
  void SampleAlongGeodesics();
  void IntegrateRadiation();
  double RadialGeodesicCoordinate(double x, double y, double z);
  void CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov);
  void ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon);
  void ContravariantGeodesicMetricDerivative(double x, double y, double z, Array<double> &dgcon);
};

#endif
