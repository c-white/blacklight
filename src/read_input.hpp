// Blacklight input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// C++ headers
#include <string>  // string

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
  ModelType model_type;
  int num_threads;

  // Data - output parameters
  OutputFormat output_format;
  std::string output_file;

  // Data - simulation parameters
  std::string simulation_file;
  double simulation_m_msun;
  double simulation_a;
  double simulation_rho_cgs;
  Coordinates simulation_coord;
  bool simulation_interp;

  // Data - plasma parameters
  double plasma_mu;
  double plasma_ne_ni;
  double plasma_rat_high;
  double plasma_rat_low;
  double plasma_sigma_max;

  // Data - formula parameters
  double formula_mass;
  double formula_spin;
  double formula_r0;
  double formula_h;
  double formula_l0;
  double formula_q;
  double formula_nup;
  double formula_cn0;
  double formula_alpha;
  double formula_a;
  double formula_beta;

  // Data - image parameters
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

  // Data - ray-tracing parameters
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

  // Functions
  void Read();
  static bool RemoveableSpace(unsigned char c);
  double ReadPole(const std::string &string, bool *p_pole_flag);
  bool ReadBool(const std::string &string);
  ModelType ReadModelType(const std::string &string);
  OutputFormat ReadOutputFormat(const std::string &string);
  Coordinates ReadCoordinates(const std::string &string);
  Camera ReadCamera(const std::string &string);
  FrequencyNormalization ReadFrequencyNormalization(const std::string &string);
  RayTerminate ReadRayTerminate(const std::string &string);
};

#endif
