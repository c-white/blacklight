// Blacklight ray tracer header

#ifndef RAY_TRACER_H_
#define RAY_TRACER_H_

// Blacklight headers
#include "array.hpp"        // Array
#include "blacklight.hpp"   // Coordinates
#include "read_athena.hpp"  // AthenaReader
#include "read_input.hpp"   // InputReader

//--------------------------------------------------------------------------------------------------

// Ray tracer
// TODO: decide if nan should be used for fallback
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
  Coordinates coord;

  // Input data - units
  double m_msun;
  double rho_unit;

  // Input data - plasma
  double plasma_mu;
  double plasma_ne_ni;
  double plasma_rat_high;
  double plasma_rat_low;
  double plasma_sigma_max;

  // Input data - image
  double im_r;
  double im_th;
  double im_ph;
  double im_rot;
  double im_width;
  int im_res;
  double im_freq;
  bool im_pole;

  // Input data - rays
  double ray_step;
  int ray_max_steps;
  bool ray_sample_interp;
  bool ray_flat;

  // Grid data
  double x1_min, x1_max, x2_min, x2_max, x3_min, x3_max;
  Array<float> x1f, x2f, x3f;
  Array<float> x1v, x2v, x3v;
  Array<float> grid_rho, grid_pgas;
  Array<float> grid_uu1, grid_uu2, grid_uu3;
  Array<float> grid_bb1, grid_bb2, grid_bb3;

  // Fallback data
  const float rho_fallback = 1.0e-6f;
  const float pgas_fallback = 1.0e-8f;
  const float uu1_fallback = 0.0f;
  const float uu2_fallback = 0.0f;
  const float uu3_fallback = 0.0f;
  const float bb1_fallback = 0.0f;
  const float bb2_fallback = 0.0f;
  const float bb3_fallback = 0.0f;

  // Sample and image data
  double r_photon;
  int im_steps;
  Array<double> im_pos, im_dir;
  Array<double> geodesic_pos, geodesic_dir, geodesic_len;
  Array<bool> sample_flags;
  Array<int> sample_num;
  Array<double> sample_pos, sample_dir, sample_len;
  Array<float> sample_rho, sample_pgas;
  Array<float> sample_uu1, sample_uu2, sample_uu3;
  Array<float> sample_bb1, sample_bb2, sample_bb3;
  Array<double> image;

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
  void CovariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcov);
  void ContravariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcon);
};

#endif
