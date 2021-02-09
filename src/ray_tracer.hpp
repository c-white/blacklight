// Blacklight ray tracer header

#ifndef RAY_TRACER_H_
#define RAY_TRACER_H_

// Blacklight headers
#include "array.hpp"        // Array
#include "blacklight.hpp"   // enumerations
#include "read_athena.hpp"  // AthenaReader
#include "read_input.hpp"   // InputReader

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

  // Input data - simulation parameters
  double simulation_m_msun;
  double simulation_rho_cgs;
  Coordinates simulation_coord;
  bool simulation_interp;
  bool simulation_block_interp;

  // Input data - plasma parameters
  double plasma_mu;
  double plasma_ne_ni;
  double plasma_rat_high;
  double plasma_rat_low;
  double plasma_sigma_max;

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

  // Input data - fallback parameters
  bool fallback_nan;
  float fallback_rho;
  float fallback_pgas;

  // Input data - image parameters
  Camera im_camera;
  double im_r;
  double im_th;
  double im_ph;
  double im_urn;
  double im_uthn;
  double im_uphn;
  double im_k_r;
  double im_k_th;
  double im_k_ph;
  double im_rot;
  double im_width;
  int im_res;
  double im_freq;
  FrequencyNormalization im_norm;
  bool im_pole;

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
  Array<float> grid_rho, grid_pgas;
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

  // External function
  void MakeImage();

  // Internal functions
  void InitializeCamera();
  void InitializeGeodesics();
  void IntegrateGeodesics();
  void TransformGeodesics();
  void SampleSimulationAlongGeodesics();
  void IntegrateSimulationRadiation();
  void IntegrateFormulaRadiation();
  double RadialGeodesicCoordinate(double x, double y, double z);
  void CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov);
  void ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon);
  void ContravariantGeodesicMetricDerivative(double x, double y, double z, Array<double> &dgcon);
  void CovariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcov);
  void ContravariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcon);
  void GeodesicSubstep(double y[9], double k[9], Array<double> &gcov, Array<double> &gcon,
      Array<double> &dgcon);
  void FindNearbyVals(int b, int k, int j, int i, double vals[8]);
};

#endif
