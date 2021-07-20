// Blacklight geodesic integrator header

#ifndef GEODESIC_INTEGRATOR_H_
#define GEODESIC_INTEGRATOR_H_

// Blacklight headers
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../utils/array.hpp"                // Array

//--------------------------------------------------------------------------------------------------

// Geodesic integrator
struct GeodesicIntegrator
{
  // Constructors and destructor
  GeodesicIntegrator(const InputReader *p_input_reader);
  GeodesicIntegrator(const GeodesicIntegrator &source) = delete;
  GeodesicIntegrator &operator=(const GeodesicIntegrator &source) = delete;
  ~GeodesicIntegrator() {}

  // Input data - general
  ModelType model_type;

  // Input data - image parameters
  Camera image_camera;
  double image_r;
  double image_th;
  double image_ph;
  double image_urn;
  double image_uthn;
  double image_uphn;
  double image_k_r;
  double image_k_th;
  double image_k_ph;
  double image_rotation;
  double image_width;
  int image_resolution;
  double image_frequency;
  FrequencyNormalization image_normalization;
  bool image_pole;

  // Input data - ray-tracing parameters
  bool ray_flat;
  RayTerminate ray_terminate;
  double ray_factor;
  double ray_step;
  int ray_max_steps;
  int ray_max_retries;
  double ray_tol_abs;
  double ray_tol_rel;
  double ray_err_factor;
  double ray_min_factor;
  double ray_max_factor;

  // Geometry data
  double bh_m;
  double bh_a;
  double r_terminate;

  // Camera data
  double momentum_factor;
  double camera_ucon[4], camera_ucov[4], camera_up_con_c[4];
  int camera_num_pix;
  Array<double> camera_pos, camera_dir;

  // Geodesic data
  int geodesic_num_steps;
  Array<double> geodesic_pos, geodesic_dir, geodesic_len;
  Array<bool> sample_flags;
  Array<int> sample_num;
  Array<double> sample_pos, sample_dir, sample_len;

  // External function
  void Integrate();

  // Internal functions - camera.cpp
  void InitializeCamera();

  // Internal functions - geodesics.cpp
  void InitializeGeodesics();
  void IntegrateGeodesics();
  void ReverseGeodesics();
  void GeodesicSubstep(double y[9], double k[9], Array<double> &gcov, Array<double> &gcon,
      Array<double> &dgcon);

  // Internal functions - geodesic_geometry.cpp
  double RadialGeodesicCoordinate(double x, double y, double z);
  void CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov);
  void ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon);
  void ContravariantGeodesicMetricDerivative(double x, double y, double z, Array<double> &dgcon);
};

#endif
