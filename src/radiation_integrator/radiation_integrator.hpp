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
  ~RadiationIntegrator() {}

  // Pointers to other objects
  const AthenaReader *p_athena_reader;

  // Input data - general
  ModelType model_type;

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
  double image_frequency;
  bool image_polarization;

  // Input data - ray-tracing parameters
  bool ray_flat;

  // Flag
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
  double camera_ucon[4], camera_ucov[4], camera_up_con_c[4];
  int camera_num_pix;
  Array<double> camera_pos, camera_dir;

  // Geodesic data
  int geodesic_num_steps;
  Array<bool> sample_flags;
  Array<int> sample_num;
  Array<double> sample_pos, sample_dir, sample_len;

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

  // Sample data
  Array<int> sample_inds;
  Array<double> sample_fracs;
  Array<bool> sample_nan, sample_fallback;
  Array<float> sample_rho, sample_pgas, sample_kappa;
  Array<float> sample_uu1, sample_uu2, sample_uu3;
  Array<float> sample_bb1, sample_bb2, sample_bb3;

  // Coefficient data
  Array<double> j_i, j_q, j_v;
  Array<double> alpha_i, alpha_q, alpha_v;
  Array<double> rho_q, rho_v;

  // Image data
  Array<double> image;

  // External function
  void Integrate(double *p_time_sample, double *p_time_integrate);

  // Internal functions - simulation_sampling.cpp
  void CalculateSimulationSampling();
  void SampleSimulation();
  void FindNearbyInds(int b, int k, int j, int i, int k_c, int j_c, int i_c, double x3, double x2,
      double x1, int inds[4]);
  float InterpolateSimple(const Array<float> &grid_vals, int b, int k, int j, int i, double f_k,
      double f_j, double f_i);
  float InterpolateAdvanced(const Array<float> &grid_vals, int m, int n);

  // Internal functions - simulation_coefficients.cpp
  void CalculateSimulationCoefficients();

  // Internal functions - formula_coefficients.cpp
  void CalculateFormulaCoefficients();

  // Internal functions - unpolarized.cpp
  void IntegrateUnpolarizedRadiation();

  // Internal functions - polarized.cpp
  void IntegratePolarizedRadiation();

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
