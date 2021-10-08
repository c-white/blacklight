// Blacklight input reader header

#ifndef INPUT_READER_H_
#define INPUT_READER_H_

// C++ headers
#include <optional>  // optional
#include <string>    // string

// Blacklight headers
#include "../blacklight.hpp"  // enums

//--------------------------------------------------------------------------------------------------

// Input reader
struct InputReader
{
  // Constructors and destructor
  InputReader(const std::string input_file_);
  InputReader(const InputReader &source) = delete;
  InputReader &operator=(const InputReader &source) = delete;
  ~InputReader();

  // Input file
  const std::string input_file;

  // Data - general
  std::optional<ModelType> model_type;
  std::optional<int> num_threads;

  // Data - output parameters
  std::optional<OutputFormat> output_format;
  std::optional<std::string> output_file;
  std::optional<bool> output_camera;

  // Data - checkpoint parameters
  std::optional<bool> checkpoint_geodesic_save;
  std::optional<bool> checkpoint_geodesic_load;
  std::optional<std::string> checkpoint_geodesic_file;
  std::optional<bool> checkpoint_sample_save;
  std::optional<bool> checkpoint_sample_load;
  std::optional<std::string> checkpoint_sample_file;

  // Data - formula parameters
  std::optional<double> formula_mass;
  std::optional<double> formula_spin;
  std::optional<double> formula_r0;
  std::optional<double> formula_h;
  std::optional<double> formula_l0;
  std::optional<double> formula_q;
  std::optional<double> formula_nup;
  std::optional<double> formula_cn0;
  std::optional<double> formula_alpha;
  std::optional<double> formula_a;
  std::optional<double> formula_beta;

  // Data - simulation parameters
  std::optional<std::string> simulation_file;
  std::optional<bool> simulation_multiple;
  std::optional<int> simulation_start;
  std::optional<int> simulation_end;
  std::optional<Coordinates> simulation_coord;
  std::optional<double> simulation_m_msun;
  std::optional<double> simulation_a;
  std::optional<double> simulation_rho_cgs;
  std::optional<std::string> simulation_kappa_name;
  std::optional<bool> simulation_interp;
  std::optional<bool> simulation_block_interp;

  // Data - plasma parameters
  std::optional<double> plasma_mu;
  std::optional<double> plasma_ne_ni;
  std::optional<PlasmaModel> plasma_model;
  std::optional<double> plasma_rat_low;
  std::optional<double> plasma_rat_high;
  std::optional<double> plasma_power_frac;
  std::optional<double> plasma_p;
  std::optional<double> plasma_gamma_min;
  std::optional<double> plasma_gamma_max;
  std::optional<double> plasma_kappa_frac;
  std::optional<double> plasma_kappa;
  std::optional<double> plasma_w;
  std::optional<double> plasma_sigma_max;

  // Data - slow light parameters
  std::optional<bool> slow_light_on;
  std::optional<bool> slow_interp;
  std::optional<int> slow_chunk_size;
  std::optional<double> slow_t_start;
  std::optional<double> slow_dt;
  std::optional<int> slow_num_images;
  std::optional<int> slow_offset;

  // Data - fallback parameters
  std::optional<bool> fallback_nan;
  std::optional<float> fallback_rho;
  std::optional<float> fallback_pgas;
  std::optional<float> fallback_kappa;

  // Data - camera parameters
  std::optional<Camera> camera_type;
  std::optional<double> camera_r;
  std::optional<double> camera_th;
  std::optional<double> camera_ph;
  std::optional<double> camera_urn;
  std::optional<double> camera_uthn;
  std::optional<double> camera_uphn;
  std::optional<double> camera_k_r;
  std::optional<double> camera_k_th;
  std::optional<double> camera_k_ph;
  std::optional<double> camera_rotation;
  std::optional<double> camera_width;
  std::optional<int> camera_resolution;
  std::optional<bool> camera_pole;

  // Data - image parameters
  std::optional<bool> image_light;
  std::optional<double> image_frequency;
  std::optional<FrequencyNormalization> image_normalization;
  std::optional<bool> image_polarization;
  std::optional<bool> image_time;
  std::optional<bool> image_length;
  std::optional<bool> image_lambda;
  std::optional<bool> image_emission;
  std::optional<bool> image_tau;
  std::optional<bool> image_lambda_ave;
  std::optional<bool> image_emission_ave;
  std::optional<bool> image_tau_int;

  // Data - rendering parameters
  std::optional<int> render_num_images;
  std::optional<int> *render_num_features = nullptr;
  std::optional<int> **render_quantities = nullptr;
  std::optional<RenderType> **render_types = nullptr;
  std::optional<double> **render_thresh_vals = nullptr;
  std::optional<double> **render_min_vals = nullptr;
  std::optional<double> **render_max_vals = nullptr;
  std::optional<double> **render_opacities = nullptr;
  std::optional<double> **render_tau_scales = nullptr;
  std::optional<double> **render_x_vals = nullptr;
  std::optional<double> **render_y_vals = nullptr;
  std::optional<double> **render_z_vals = nullptr;

  // Data - ray-tracing parameters
  std::optional<bool> ray_flat;
  std::optional<RayTerminate> ray_terminate;
  std::optional<double> ray_factor;
  std::optional<double> ray_step;
  std::optional<int> ray_max_steps;
  std::optional<int> ray_max_retries;
  std::optional<double> ray_tol_abs;
  std::optional<double> ray_tol_rel;
  std::optional<double> ray_err_factor;
  std::optional<double> ray_min_factor;
  std::optional<double> ray_max_factor;

  // Data - adaptive parameters
  std::optional<bool> adaptive_on;
  std::optional<int> adaptive_block_size;
  std::optional<int> adaptive_max_level;
  std::optional<double> adaptive_val_cut;
  std::optional<double> adaptive_val_frac;
  std::optional<double> adaptive_abs_grad_cut;
  std::optional<double> adaptive_abs_grad_frac;
  std::optional<double> adaptive_rel_grad_cut;
  std::optional<double> adaptive_rel_grad_frac;
  std::optional<double> adaptive_abs_lapl_cut;
  std::optional<double> adaptive_abs_lapl_frac;
  std::optional<double> adaptive_rel_lapl_cut;
  std::optional<double> adaptive_rel_lapl_frac;

  // External functions
  int Read();

  // Internal functions - input_reader.cpp
  static bool RemoveableSpace(unsigned char c);
  double ReadPole(const std::string &string, std::optional<bool> *p_pole_flag);
  bool ReadBool(const std::string &string);
  void ReadRender(const std::string &key, const std::string &val);

  // Internal functions - enum_readers.cpp
  ModelType ReadModelType(const std::string &string);
  OutputFormat ReadOutputFormat(const std::string &string);
  Coordinates ReadCoordinates(const std::string &string);
  PlasmaModel ReadPlasmaModel(const std::string &string);
  Camera ReadCamera(const std::string &string);
  FrequencyNormalization ReadFrequencyNormalization(const std::string &string);
  RayTerminate ReadRayTerminate(const std::string &string);
};

#endif
