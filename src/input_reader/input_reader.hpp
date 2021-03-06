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
  ~InputReader() {}

  // Input file
  const std::string input_file;

  // Data - general
  std::optional<ModelType> model_type;
  std::optional<int> num_threads;

  // Data - output parameters
  std::optional<OutputFormat> output_format;
  std::optional<std::string> output_file_template;
  std::optional<bool> output_params;
  std::optional<bool> output_camera;

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
  std::optional<std::string> simulation_file_template;
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
  std::optional<double> plasma_rat_high;
  std::optional<double> plasma_rat_low;
  std::optional<double> plasma_sigma_max;

  // Data - fallback parameters
  std::optional<bool> fallback_nan;
  std::optional<float> fallback_rho;
  std::optional<float> fallback_pgas;
  std::optional<float> fallback_kappa;

  // Data - image parameters
  std::optional<Camera> image_camera;
  std::optional<double> image_r;
  std::optional<double> image_th;
  std::optional<double> image_ph;
  std::optional<double> image_urn;
  std::optional<double> image_uthn;
  std::optional<double> image_uphn;
  std::optional<double> image_k_r;
  std::optional<double> image_k_th;
  std::optional<double> image_k_ph;
  std::optional<double> image_rotation;
  std::optional<double> image_width;
  std::optional<int> image_resolution;
  std::optional<double> image_frequency;
  std::optional<FrequencyNormalization> image_normalization;
  std::optional<bool> image_polarization;
  std::optional<bool> image_pole;

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

  // String substitution data
  bool multiple_runs;
  int num_runs;
  std::string simulation_file_formatted;
  std::string::size_type simulation_pos_open;
  std::string::size_type simulation_pos_close;
  int simulation_field_length;
  std::string output_file_formatted;
  std::string::size_type output_pos_open;
  std::string::size_type output_pos_close;
  int output_field_length;

  // External functions
  int Read();
  void AdjustFileNames(int index);

  // Internal functions
  static bool RemoveableSpace(unsigned char c);
  double ReadPole(const std::string &string, std::optional<bool> *p_pole_flag);
  bool ReadBool(const std::string &string);
  ModelType ReadModelType(const std::string &string);
  OutputFormat ReadOutputFormat(const std::string &string);
  Coordinates ReadCoordinates(const std::string &string);
  PlasmaModel ReadPlasmaModel(const std::string &string);
  Camera ReadCamera(const std::string &string);
  FrequencyNormalization ReadFrequencyNormalization(const std::string &string);
  RayTerminate ReadRayTerminate(const std::string &string);
};

#endif
