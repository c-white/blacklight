// Blacklight input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// C++ headers
#include <optional>  // optional
#include <string>    // string

// Blacklight headers
#include "blacklight.hpp"  // enumerations

//--------------------------------------------------------------------------------------------------

// Input reader
struct InputReader
{
  // Constructor
  InputReader(const std::string input_file_);

  // Input file
  const std::string input_file;

  // Data - general
  std::optional<ModelType> model_type;
  std::optional<int> num_threads;

  // Data - output parameters
  std::optional<OutputFormat> output_format;
  std::optional<std::string> output_file;

  // Data - simulation parameters
  std::optional<std::string> simulation_file;
  std::optional<double> simulation_m_msun;
  std::optional<double> simulation_a;
  std::optional<double> simulation_rho_cgs;
  std::optional<Coordinates> simulation_coord;
  std::optional<bool> simulation_interp;
  std::optional<bool> simulation_block_interp;

  // Data - plasma parameters
  std::optional<double> plasma_mu;
  std::optional<double> plasma_ne_ni;
  std::optional<double> plasma_rat_high;
  std::optional<double> plasma_rat_low;
  std::optional<double> plasma_sigma_max;

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

  // Data - fallback parameters
  std::optional<bool> fallback_nan;
  std::optional<float> fallback_rho;
  std::optional<float> fallback_pgas;

  // Data - image parameters
  std::optional<Camera> im_camera;
  std::optional<double> im_r;
  std::optional<double> im_th;
  std::optional<double> im_ph;
  std::optional<double> im_urn;
  std::optional<double> im_uthn;
  std::optional<double> im_uphn;
  std::optional<double> im_k_r;
  std::optional<double> im_k_th;
  std::optional<double> im_k_ph;
  std::optional<double> im_rot;
  std::optional<double> im_width;
  std::optional<int> im_res;
  std::optional<double> im_freq;
  std::optional<FrequencyNormalization> im_norm;
  std::optional<bool> im_pole;

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

  // Functions
  void Read();
  static bool RemoveableSpace(unsigned char c);
  double ReadPole(const std::string &string, std::optional<bool> *p_pole_flag);
  bool ReadBool(const std::string &string);
  ModelType ReadModelType(const std::string &string);
  OutputFormat ReadOutputFormat(const std::string &string);
  Coordinates ReadCoordinates(const std::string &string);
  Camera ReadCamera(const std::string &string);
  FrequencyNormalization ReadFrequencyNormalization(const std::string &string);
  RayTerminate ReadRayTerminate(const std::string &string);
};

#endif
