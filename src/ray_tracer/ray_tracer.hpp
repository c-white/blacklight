// Blacklight ray tracer header

#ifndef RAY_TRACER_H_
#define RAY_TRACER_H_

// Blacklight headers
#include "../blacklight.hpp"                   // enums
#include "../athena_reader/athena_reader.hpp"  // AthenaReader
#include "../input_reader/input_reader.hpp"    // InputReader
#include "../utils/array.hpp"                  // Array

//--------------------------------------------------------------------------------------------------

// Ray tracer
struct RayTracer
{
  // Constructors and destructor
  RayTracer(const InputReader *p_input_reader, const AthenaReader *p_athena_reader);
  RayTracer(const RayTracer &source) = delete;
  RayTracer &operator=(const RayTracer &source) = delete;
  ~RayTracer() {}

  // Input data - general
  ModelType model_type;

  // Parameters
  double bh_m;
  double bh_a;

  // Input data - formula parameters
  double formula_mass;
  double formula_r0;
  double formula_h;
  double formula_l0;
  double formula_q;
  double formula_nup;
  double formula_cn0;
  double formula_alpha;
  double formula_a;
  double formula_beta;

  // Input data - simulation parameters
  Coordinates simulation_coord;
  double simulation_m_msun;
  double simulation_rho_cgs;
  bool simulation_interp;
  bool simulation_block_interp;

  // Input data - plasma parameters
  double plasma_mu;
  double plasma_ne_ni;
  PlasmaModel plasma_model;
  double plasma_rat_high;
  double plasma_rat_low;
  double plasma_sigma_max;

  // Input data - fallback parameters
  bool fallback_nan;
  float fallback_rho;
  float fallback_pgas;
  float fallback_kappa;

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
  bool image_polarization;
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

  // Grid data
  int n_3_root;
  int max_level;
  Array<int> n_3_level;
  Array<int> levels, locations;
  Array<float> x1f, x2f, x3f;
  Array<float> x1v, x2v, x3v;
  Array<float> grid_rho, grid_pgas, grid_kappa;
  Array<float> grid_uu1, grid_uu2, grid_uu3;
  Array<float> grid_bb1, grid_bb2, grid_bb3;

  // Fallback data
  const float fallback_uu1 = 0.0f;
  const float fallback_uu2 = 0.0f;
  const float fallback_uu3 = 0.0f;
  const float fallback_bb1 = 0.0f;
  const float fallback_bb2 = 0.0f;
  const float fallback_bb3 = 0.0f;

  // Limiter data
  const double delta_tau_max = 100.0;

  // Sample and image data
  double mass_msun;
  double r_terminate;
  double momentum_factor;
  int image_steps;
  double camera_ucon[4], camera_ucov[4], camera_up_con_c[4];
  Array<double> image_position, image_direction;
  Array<double> geodesic_pos, geodesic_dir, geodesic_len;
  Array<bool> sample_flags;
  Array<int> sample_num;
  Array<double> sample_pos, sample_dir, sample_len;
  Array<float> sample_rho, sample_pgas, sample_kappa;
  Array<float> sample_uu1, sample_uu2, sample_uu3;
  Array<float> sample_bb1, sample_bb2, sample_bb3;
  Array<double> image;

  // External function
  void MakeImage();

  // Internal functions
  void InitializeCamera();
  void InitializeGeodesics();
  void IntegrateGeodesics();
  void TransformGeodesics();
  void GeodesicSubstep(double y[9], double k[9], Array<double> &gcov, Array<double> &gcon,
      Array<double> &dgcon);
  void SampleSimulationAlongGeodesics();
  void IntegrateSimulationUnpolarizedRadiation();
  void IntegrateSimulationPolarizedRadiation();
  void FindNearbyVals(int b, int k, int j, int i, double vals[8]);
  void IntegrateFormulaRadiation();
  double RadialGeodesicCoordinate(double x, double y, double z);
  void CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov);
  void ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon);
  void ContravariantGeodesicMetricDerivative(double x, double y, double z, Array<double> &dgcon);
  void CovariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcov);
  void ContravariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcon);
  void CoordinateConnection(double x1, double x2, double x3, Array<double> &connection);
  void GeodesicCoordinateJacobian(double x, double y, double z, Array<double> &jacobian);
  void Tetrad(const double ucon[4], const double ucov[4], const double kcon[4],
      const double kcov[4], const double up_con[4], const Array<double> &gcov,
      const Array<double> &gcon, Array<double> &tetrad);
};

#endif
