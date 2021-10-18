// Blacklight radiation integrator header

#ifndef RADIATION_INTEGRATOR_H_
#define RADIATION_INTEGRATOR_H_

// Blacklight headers
#include "../blacklight.hpp"                               // enums
#include "../athena_reader/athena_reader.hpp"              // AthenaReader
#include "../geodesic_integrator/geodesic_integrator.hpp"  // GeodesicIntegrator
#include "../input_reader/input_reader.hpp"                // InputReader
#include "../utils/array.hpp"                              // Array

//--------------------------------------------------------------------------------------------------

// Radiation integrator
struct RadiationIntegrator
{
  // Constructors and destructor
  RadiationIntegrator(const InputReader *p_input_reader,
      const GeodesicIntegrator *p_geodesic_integrator, const AthenaReader *p_athena_reader_);
  RadiationIntegrator(const RadiationIntegrator &source) = delete;
  RadiationIntegrator &operator=(const RadiationIntegrator &source) = delete;
  ~RadiationIntegrator();

  // Pointers to other objects
  const AthenaReader *p_athena_reader;

  // Input data - general
  ModelType model_type;
  int num_threads;

  // Input data - checkpoints
  bool checkpoint_sample_save;
  bool checkpoint_sample_load;
  std::string checkpoint_sample_file;

  // Input data - simulation parameters
  Coordinates simulation_coord;
  double simulation_m_msun;
  double simulation_rho_cgs;
  bool simulation_interp;
  bool simulation_block_interp;

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

  // Input data - camera parameters
  double camera_r;
  int camera_resolution;

  // Input data - ray-tracing parameters
  bool ray_flat;

  // Input data - image parameters
  bool image_light;
  double image_frequency;
  bool image_polarization;
  bool image_time;
  bool image_length;
  bool image_lambda;
  bool image_emission;
  bool image_tau;
  bool image_lambda_ave;
  bool image_emission_ave;
  bool image_tau_int;

  // Input data - rendering parameters
  int render_num_images;
  int *render_num_features = nullptr;
  int **render_quantities = nullptr;
  RenderType **render_types = nullptr;
  double **render_min_vals = nullptr;
  double **render_max_vals = nullptr;
  double **render_thresh_vals = nullptr;
  double **render_tau_scales = nullptr;
  double **render_opacities = nullptr;
  double **render_x_vals = nullptr;
  double **render_y_vals = nullptr;
  double **render_z_vals = nullptr;

  // Input data - slow-light parameters
  bool slow_light_on;
  bool slow_interp;
  int slow_chunk_size;
  double slow_t_start;
  double slow_dt;

  // Input data - adaptive parameters
  int adaptive_max_level;
  int adaptive_block_size;
  double adaptive_val_cut;
  double adaptive_val_frac;
  double adaptive_abs_grad_cut;
  double adaptive_abs_grad_frac;
  double adaptive_rel_grad_cut;
  double adaptive_rel_grad_frac;
  double adaptive_abs_lapl_cut;
  double adaptive_abs_lapl_frac;
  double adaptive_rel_lapl_cut;
  double adaptive_rel_lapl_frac;

  // Input data - plasma parameters
  double plasma_mu;
  double plasma_ne_ni;
  double plasma_thermal_frac;
  PlasmaModel plasma_model;
  double plasma_rat_low;
  double plasma_rat_high;
  double plasma_power_frac;
  double plasma_p;
  double plasma_gamma_min;
  double plasma_gamma_max;
  double plasma_kappa_frac;
  double plasma_kappa;
  double plasma_w;

  // Input data - cut parameters
  double cut_rho_min;
  double cut_rho_max;
  double cut_n_e_min;
  double cut_n_e_max;
  double cut_p_gas_min;
  double cut_p_gas_max;
  double cut_theta_e_min;
  double cut_theta_e_max;
  double cut_b_min;
  double cut_b_max;
  double cut_sigma_min;
  double cut_sigma_max;
  double cut_beta_inverse_min;
  double cut_beta_inverse_max;
  bool cut_omit_near;
  bool cut_omit_far;
  double cut_omit_in;
  double cut_omit_out;
  bool cut_plane;
  double cut_plane_origin_x, cut_plane_origin_y, cut_plane_origin_z;
  double cut_plane_normal_x, cut_plane_normal_y, cut_plane_normal_z;

  // Input data - fallback parameters
  bool fallback_nan;
  float fallback_rho;
  float fallback_pgas;
  float fallback_kappa;

  // Flag for tracking function calls
  bool first_time = true;

  // Geometry data
  double bh_m;
  double bh_a;
  double mass_msun;

  // Fallback data
  const float fallback_uu1 = 0.0f;
  const float fallback_uu2 = 0.0f;
  const float fallback_uu3 = 0.0f;
  const float fallback_bb1 = 0.0f;
  const float fallback_bb2 = 0.0f;
  const float fallback_bb3 = 0.0f;

  // Limiter data
  const double delta_tau_max = 100.0;

  // Camera data
  double momentum_factor;
  double camera_x[4];
  double camera_u_con[4], camera_u_cov[4];
  double camera_vert_con_c[4];
  int camera_num_pix;
  Array<double> *camera_pos = nullptr;
  Array<double> *camera_dir = nullptr;

  // Geodesic data
  int *geodesic_num_steps = nullptr;
  Array<bool> *sample_flags = nullptr;
  Array<int> *sample_num = nullptr;
  Array<double> *sample_pos = nullptr;
  Array<double> *sample_dir = nullptr;
  Array<double> *sample_len = nullptr;

  // Grid data
  int n_3_root;
  int max_level;
  Array<int> n_3_level;
  Array<int> levels;
  Array<int> locations;
  Array<float> x1f, x2f, x3f;
  Array<float> x1v, x2v, x3v;
  float *time;
  Array<float> *grid_prim;
  Array<float> *grid_bb;
  int ind_rho, ind_pgas, ind_kappa;
  int ind_uu1, ind_uu2, ind_uu3;
  int ind_bb1, ind_bb2, ind_bb3;

  // Sample data
  Array<int> *sample_inds = nullptr;
  Array<double> *sample_fracs = nullptr;
  Array<bool> *sample_nan = nullptr;
  Array<bool> *sample_cut = nullptr;
  Array<bool> *sample_fallback = nullptr;
  Array<float> *sample_rho = nullptr;
  Array<float> *sample_pgas = nullptr;
  Array<float> *sample_kappa = nullptr;
  Array<float> *sample_uu1 = nullptr;
  Array<float> *sample_uu2 = nullptr;
  Array<float> *sample_uu3 = nullptr;
  Array<float> *sample_bb1 = nullptr;
  Array<float> *sample_bb2 = nullptr;
  Array<float> *sample_bb3 = nullptr;
  float extrapolation_tolerance;

  // Coefficient data
  Array<double> *j_i = nullptr;
  Array<double> *j_q = nullptr;
  Array<double> *j_v = nullptr;
  Array<double> *alpha_i = nullptr;
  Array<double> *alpha_q = nullptr;
  Array<double> *alpha_v = nullptr;
  Array<double> *rho_q = nullptr;
  Array<double> *rho_v = nullptr;
  Array<double> *cell_values = nullptr;

  // Image data
  Array<double> *image = nullptr;
  int image_num_quantities = 0;
  int image_offset_time = 0;
  int image_offset_length = 0;
  int image_offset_lambda = 0;
  int image_offset_emission = 0;
  int image_offset_tau = 0;
  int image_offset_lambda_ave = 0;
  int image_offset_emission_ave = 0;
  int image_offset_tau_int = 0;

  // Rendering data
  Array<double> *render = nullptr;

  // Adaptive data
  int adaptive_level = 0;
  int adaptive_num_levels;
  int linear_root_blocks;
  int block_num_pix;
  int *block_counts = nullptr;
  Array<bool> *refinement_flags = nullptr;
  Array<double> *image_blocks = nullptr;

  // Precalculated values
  double power_jj, power_jj_q, power_jj_v;
  double power_aa, power_aa_q, power_aa_v;
  double power_rho, power_rho_q, power_rho_v;
  double kappa_jj_low, kappa_jj_low_q, kappa_jj_low_v;
  double kappa_jj_high, kappa_jj_high_q, kappa_jj_high_v;
  double kappa_jj_x_i, kappa_jj_x_q, kappa_jj_x_v;
  double kappa_aa_low, kappa_aa_low_q, kappa_aa_low_v;
  double kappa_aa_high, kappa_aa_high_i, kappa_aa_high_q, kappa_aa_high_v;
  double kappa_aa_x_i, kappa_aa_x_q, kappa_aa_x_v;
  double kappa_rho_frac;
  double kappa_rho_q_low_a, kappa_rho_q_low_b, kappa_rho_q_low_c, kappa_rho_q_low_d,
      kappa_rho_q_low_e;
  double kappa_rho_q_high_a, kappa_rho_q_high_b, kappa_rho_q_high_c, kappa_rho_q_high_d,
      kappa_rho_q_high_e;
  double kappa_rho_v;
  double kappa_rho_v_low_a, kappa_rho_v_low_b;
  double kappa_rho_v_high_a, kappa_rho_v_high_b;

  // External function
  bool Integrate(int snapshot, double *p_time_sample, double *p_time_image, double *p_time_render);

  // Internal functions - sample_checkpoint.cpp
  void SaveSampling();
  void LoadSampling();

  // Internal functions - simulation_sampling.cpp
  void ObtainGridData();
  void CalculateSimulationSampling(int snapshot);
  void SampleSimulation();
  void FindNearbyInds(int b, int k, int j, int i, int k_c, int j_c, int i_c, double x3, double x2,
      double x1, int inds[4]);
  double InterpolateSimple(const Array<float> &grid_vals, int grid_ind, int b, int k, int j, int i,
      double f_k, double f_j, double f_i);
  double InterpolateAdvanced(const Array<float> &grid_vals, int grid_ind, int m, int n);

  // Internal functions - simulation_coefficients.cpp
  void CalculateSimulationCoefficients();
  double Hypergeometric(double alpha, double beta, double gamma, double z);

  // Internal functions - formula_coefficients.cpp
  void CalculateFormulaCoefficients();

  // Internal functions - unpolarized.cpp
  void IntegrateUnpolarizedRadiation();

  // Internal functions - polarized.cpp
  void IntegratePolarizedRadiation();

  // Internal functions - rendering.cpp
  void Render();

  // Internal functions - radiation_adaptive.cpp
  bool CheckAdaptiveRefinement();
  bool EvaluateBlock(int thread);

  // Internal functions - radiation_geometry.cpp
  double RadialGeodesicCoordinate(double x, double y, double z);
  void CKSToSKS(double *p_x1, double *p_x2, double *p_x3);
  void CoordinateJacobian(double x, double y, double z, Array<double> &jacobian);
  void CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov);
  void ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon);
  void GeodesicConnection(double x, double y, double z, Array<double> &connection);
  void CovariantSimulationMetric(double x, double y, double z, Array<double> &gcov);
  void ContravariantSimulationMetric(double x, double y, double z, Array<double> &gcon);
  void Tetrad(const double ucon[4], const double ucov[4], const double kcon[4],
      const double kcov[4], const double up_con[4], const Array<double> &gcov,
      const Array<double> &gcon, Array<double> &tetrad);
};

#endif
