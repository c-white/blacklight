// Blacklight input reader

// C++ headers
#include <algorithm>  // remove_if
#include <cctype>     // isspace
#include <fstream>    // ifstream
#include <optional>   // optional
#include <sstream>    // ostringstream
#include <string>     // getline, stod, stof, stoi, string

// Blacklight headers
#include "input_reader.hpp"
#include "../blacklight.hpp"        // math, enums
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Input reader constructor
// Inputs:
//   input_file_: name of input file
InputReader::InputReader(const std::string input_file_)
  : input_file(input_file_) {}

//--------------------------------------------------------------------------------------------------

// Input reader read and initialize function
// Inputs: (none)
// Outputs:
//   returned value: number of runs to perform based on inputs
int InputReader::Read()
{
  // Open input file
  std::ifstream input_stream(input_file);
  if (not input_stream.is_open())
    throw BlacklightException("Could not open input file.");

  // Process file line by line
  for (std::string line; std::getline(input_stream, line); )
  {
    // Remove spaces
    line.erase(std::remove_if(line.begin(), line.end(), RemoveableSpace), line.end());

    // Remove comments
    std::string::size_type pos = line.find('#');
    if (pos != std::string::npos)
      line.erase(pos);

    // Skip blank lines
    if (line.empty())
      continue;

    // Split on '='
    pos = line.find('=');
    if (pos == std::string::npos)
      throw BlacklightException("Invalid assignment in input file.");
    std::string key = line.substr(0, pos);
    std::string val = line.substr(pos + 1, line.size());

    // Store general data
    if (key == "model_type")
      model_type = ReadModelType(val);
    else if (key == "num_threads")
      num_threads = std::stoi(val);

    // Store output parameters
    else if (key == "output_format")
      output_format = ReadOutputFormat(val);
    else if (key == "output_file")
      output_file = val;
    else if (key == "output_params")
      output_params = ReadBool(val);
    else if (key == "output_camera")
      output_camera = ReadBool(val);

    // Store checkpoint parameters
    else if (key == "checkpoint_geodesic_save")
      checkpoint_geodesic_save = ReadBool(val);
    else if (key == "checkpoint_geodesic_load")
      checkpoint_geodesic_load = ReadBool(val);
    else if (key == "checkpoint_geodesic_file")
      checkpoint_geodesic_file = val;
    else if (key == "checkpoint_sample_save")
      checkpoint_sample_save = ReadBool(val);
    else if (key == "checkpoint_sample_load")
      checkpoint_sample_load = ReadBool(val);
    else if (key == "checkpoint_sample_file")
      checkpoint_sample_file = val;

    // Store formula parameters
    else if (key == "formula_mass")
      formula_mass = std::stod(val);
    else if (key == "formula_spin")
      formula_spin = std::stod(val);
    else if (key == "formula_r0")
      formula_r0 = std::stod(val);
    else if (key == "formula_h")
      formula_h = std::stod(val);
    else if (key == "formula_l0")
      formula_l0 = std::stod(val);
    else if (key == "formula_q")
      formula_q = std::stod(val);
    else if (key == "formula_nup")
      formula_nup = std::stod(val);
    else if (key == "formula_cn0")
      formula_cn0 = std::stod(val);
    else if (key == "formula_alpha")
      formula_alpha = std::stod(val);
    else if (key == "formula_a")
      formula_a = std::stod(val);
    else if (key == "formula_beta")
      formula_beta = std::stod(val);

    // Store simulation parameters
    else if (key == "simulation_file")
      simulation_file = val;
    else if (key == "simulation_multiple")
      simulation_multiple = ReadBool(val);
    else if (key == "simulation_start")
      simulation_start = std::stoi(val);
    else if (key == "simulation_end")
      simulation_end = std::stoi(val);
    else if (key == "simulation_coord")
      simulation_coord = ReadCoordinates(val);
    else if (key == "simulation_m_msun")
      simulation_m_msun = std::stod(val);
    else if (key == "simulation_a")
      simulation_a = std::stod(val);
    else if (key == "simulation_rho_cgs")
      simulation_rho_cgs = std::stod(val);
    else if (key == "simulation_kappa_name")
      simulation_kappa_name = val;
    else if (key == "simulation_interp")
      simulation_interp = ReadBool(val);
    else if (key == "simulation_block_interp")
      simulation_block_interp = ReadBool(val);

    // Store plasma parameters
    else if (key == "plasma_mu")
      plasma_mu = std::stod(val);
    else if (key == "plasma_ne_ni")
      plasma_ne_ni = std::stod(val);
    else if (key == "plasma_model")
      plasma_model = ReadPlasmaModel(val);
    else if (key == "plasma_rat_low")
      plasma_rat_low = std::stod(val);
    else if (key == "plasma_rat_high")
      plasma_rat_high = std::stod(val);
    else if (key == "plasma_power_frac")
      plasma_power_frac = std::stod(val);
    else if (key == "plasma_p")
      plasma_p = std::stod(val);
    else if (key == "plasma_gamma_min")
      plasma_gamma_min = std::stod(val);
    else if (key == "plasma_gamma_max")
      plasma_gamma_max = std::stod(val);
    else if (key == "plasma_kappa_frac")
      plasma_kappa_frac = std::stod(val);
    else if (key == "plasma_kappa")
      plasma_kappa = std::stod(val);
    else if (key == "plasma_w")
      plasma_w = std::stod(val);
    else if (key == "plasma_sigma_max")
      plasma_sigma_max = std::stod(val);

    // Store slow light parameters
    else if (key == "slow_light_on")
      slow_light_on = ReadBool(val);
    else if (key == "slow_interp")
      slow_interp = ReadBool(val);
    else if (key == "slow_chunk_size")
      slow_chunk_size = std::stoi(val);
    else if (key == "slow_t_start")
      slow_t_start = std::stod(val);
    else if (key == "slow_dt")
      slow_dt = std::stod(val);
    else if (key == "slow_num_images")
      slow_num_images = std::stoi(val);
    else if (key == "slow_offset")
      slow_offset = std::stoi(val);

    // Store fallback parameters
    else if (key == "fallback_nan")
      fallback_nan = ReadBool(val);
    else if (key == "fallback_rho")
      fallback_rho = std::stof(val);
    else if (key == "fallback_pgas")
      fallback_pgas = std::stof(val);
    else if (key == "fallback_kappa")
      fallback_kappa = std::stof(val);

    // Store camera parameters
    else if (key == "camera_type")
      camera_type = ReadCamera(val);
    else if (key == "camera_r")
      camera_r = std::stod(val);
    else if (key == "camera_th")
      camera_th = ReadPole(val, &camera_pole) * math::pi / 180.0;
    else if (key == "camera_ph")
      camera_ph = std::stod(val) * math::pi / 180.0;
    else if (key == "camera_urn")
      camera_urn = std::stod(val);
    else if (key == "camera_uthn")
      camera_uthn = std::stod(val);
    else if (key == "camera_uphn")
      camera_uphn = std::stod(val);
    else if (key == "camera_k_r")
      camera_k_r = std::stod(val);
    else if (key == "camera_k_th")
      camera_k_th = std::stod(val);
    else if (key == "camera_k_ph")
      camera_k_ph = std::stod(val);
    else if (key == "camera_rotation")
      camera_rotation = std::stod(val) * math::pi / 180.0;
    else if (key == "camera_width")
      camera_width = std::stod(val);
    else if (key == "camera_resolution")
      camera_resolution = std::stoi(val);

    // Store image parameters
    else if (key == "image_frequency")
      image_frequency = std::stod(val);
    else if (key == "image_normalization")
      image_normalization = ReadFrequencyNormalization(val);
    else if (key == "image_polarization")
      image_polarization = ReadBool(val);

    // Store ray-tracing parameters
    else if (key == "ray_flat")
      ray_flat = ReadBool(val);
    else if (key == "ray_terminate")
      ray_terminate = ReadRayTerminate(val);
    else if (key == "ray_factor")
      ray_factor = std::stod(val);
    else if (key == "ray_step")
      ray_step = std::stod(val);
    else if (key == "ray_max_steps")
      ray_max_steps = std::stoi(val);
    else if (key == "ray_max_retries")
      ray_max_retries = std::stoi(val);
    else if (key == "ray_tol_abs")
      ray_tol_abs = std::stod(val);
    else if (key == "ray_tol_rel")
      ray_tol_rel = std::stod(val);
    else if (key == "ray_err_factor")
      ray_err_factor = std::stod(val);
    else if (key == "ray_min_factor")
      ray_min_factor = std::stod(val);
    else if (key == "ray_max_factor")
      ray_max_factor = std::stod(val);

    // Store adaptive parameters
    else if (key == "adaptive_on")
      adaptive_on = ReadBool(val);
    else if (key == "adaptive_block_size")
      adaptive_block_size = std::stoi(val);
    else if (key == "adaptive_max_level")
      adaptive_max_level = std::stoi(val);
    else if (key == "adaptive_val_cut")
      adaptive_val_cut = std::stod(val);
    else if (key == "adaptive_val_frac")
      adaptive_val_frac = std::stod(val);
    else if (key == "adaptive_abs_grad_cut")
      adaptive_abs_grad_cut = std::stod(val);
    else if (key == "adaptive_abs_grad_frac")
      adaptive_abs_grad_frac = std::stod(val);
    else if (key == "adaptive_rel_grad_cut")
      adaptive_rel_grad_cut = std::stod(val);
    else if (key == "adaptive_rel_grad_frac")
      adaptive_rel_grad_frac = std::stod(val);
    else if (key == "adaptive_abs_lapl_cut")
      adaptive_abs_lapl_cut = std::stod(val);
    else if (key == "adaptive_abs_lapl_frac")
      adaptive_abs_lapl_frac = std::stod(val);
    else if (key == "adaptive_rel_lapl_cut")
      adaptive_rel_lapl_cut = std::stod(val);
    else if (key == "adaptive_rel_lapl_frac")
      adaptive_rel_lapl_frac = std::stod(val);

    // Handle unknown entry
    else
    {
      std::ostringstream message;
      message << "Unknown key (" << val << ") in input file.";
      throw BlacklightException(message.str().c_str());
    }
  }

  // Count number of runs to do
  int num_runs = 1;
  if (model_type.value() == ModelType::simulation and simulation_multiple.value())
  {
    if (slow_light_on.value())
      num_runs = slow_num_images.value();
    else
      num_runs = simulation_end.value() - simulation_start.value() + 1;
  }
  return num_runs;
}

//--------------------------------------------------------------------------------------------------

// Definition of what constitutes a space
// Inputs:
//   c: character to be tested
// Outputs:
//   returned value: flag indicating character is a removable space
bool InputReader::RemoveableSpace(unsigned char c)
{
  return std::isspace(c) != 0;
}

//--------------------------------------------------------------------------------------------------

// Function for converting a string and determining if it exactly matches angles representing poles
// Inputs:
//   string: string to be converted
// Outputs:
//   returned value: double conversion of string
//   p_pole_flag: value set to true only if input string converts to 0.0 or 180.0
double InputReader::ReadPole(const std::string &string, std::optional<bool> *p_pole_flag)
{
  double val = std::stod(string);
  if (val == 0.0 or val == 180.0)
    *p_pole_flag = true;
  else
    *p_pole_flag = false;
  return val;
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as booleans
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: true or false
// Notes:
//   "true" evaluates to true, "false" to false, and anything else throws and exception.
bool InputReader::ReadBool(const std::string &string)
{
  if (string == "true")
    return true;
  else if (string == "false")
    return false;
  else
    throw BlacklightException("Unknown string used for boolean value.");
}
