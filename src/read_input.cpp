// Blacklight input reader

// C++ headers
#include <algorithm>  // remove_if
#include <cctype>     // isspace
#include <fstream>    // ifstream
#include <optional>   // optional
#include <string>     // getline, stod, stof, stoi, string

// Blacklight headers
#include "read_input.hpp"
#include "blacklight.hpp"  // math, enumerations
#include "exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Input reader constructor
// Inputs:
//   input_file_: name of input file
InputReader::InputReader(const std::string input_file_)
  : input_file(input_file_) {}

//--------------------------------------------------------------------------------------------------

// Input reader read and initialize function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Initializes all member objects.
void InputReader::Read()
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
    else if (key == "simulation_m_msun")
      simulation_m_msun = std::stod(val);
    else if (key == "simulation_a")
      simulation_a = std::stod(val);
    else if (key == "simulation_rho_cgs")
      simulation_rho_cgs = std::stod(val);
    else if (key == "simulation_coord")
      simulation_coord = ReadCoordinates(val);
    else if (key == "simulation_interp")
      simulation_interp = ReadBool(val);
    else if (key == "simulation_block_interp")
      simulation_block_interp = ReadBool(val);

    // Store plasma parameters
    else if (key == "plasma_mu")
      plasma_mu = std::stod(val);
    else if (key == "plasma_ne_ni")
      plasma_ne_ni = std::stod(val);
    else if (key == "plasma_rat_high")
      plasma_rat_high = std::stod(val);
    else if (key == "plasma_rat_low")
      plasma_rat_low = std::stod(val);
    else if (key == "plasma_sigma_max")
      plasma_sigma_max = std::stod(val);

    // Store fallback parameters
    else if (key == "fallback_nan")
      fallback_nan = ReadBool(val);
    else if (key == "fallback_rho")
      fallback_rho = std::stof(val);
    else if (key == "fallback_pgas")
      fallback_pgas = std::stof(val);

    // Store image parameters
    else if (key == "image_camera")
      image_camera = ReadCamera(val);
    else if (key == "image_r")
      image_r = std::stod(val);
    else if (key == "image_th")
      image_th = ReadPole(val, &image_pole) * math::pi / 180.0;
    else if (key == "image_ph")
      image_ph = std::stod(val) * math::pi / 180.0;
    else if (key == "image_urn")
      image_urn = std::stod(val);
    else if (key == "image_uthn")
      image_uthn = std::stod(val);
    else if (key == "image_uphn")
      image_uphn = std::stod(val);
    else if (key == "image_k_r")
      image_k_r = std::stod(val);
    else if (key == "image_k_th")
      image_k_th = std::stod(val);
    else if (key == "image_k_ph")
      image_k_ph = std::stod(val);
    else if (key == "image_rotation")
      image_rotation = std::stod(val) * math::pi / 180.0;
    else if (key == "image_width")
      image_width = std::stod(val);
    else if (key == "image_resolution")
      image_resolution = std::stoi(val);
    else if (key == "image_frequency")
      image_frequency = std::stod(val);
    else if (key == "image_normalization")
      image_normalization = ReadFrequencyNormalization(val);

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

    // Handle unknown entry
    else
      throw BlacklightException("Unknown key in input file.");
  }
  return;
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

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as ModelType enumerations
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid ModelType
// Notes:
//   Valid options:
//     "simulation": Athena++ output
//     "formula": parameterized formula from 2020 ApJ 897 148
ModelType InputReader::ReadModelType(const std::string &string)
{
  if (string == "simulation")
    return ModelType::simulation;
  else if (string == "formula")
    return ModelType::formula;
  else
    throw BlacklightException("Unknown string used for ModelType value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as OutputFormat enumerations
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid OutputFormat
// Notes:
//   Valid options:
//     "raw": raw binary values of image array with no extra data
//     "npy": NumPy .npy file with image array and minimal metadata
//     "npz": NumPy .npz file with image and metadata arrays
OutputFormat InputReader::ReadOutputFormat(const std::string &string)
{
  if (string == "raw")
    return OutputFormat::raw;
  else if (string == "npy")
    return OutputFormat::npy;
  else if (string == "npz")
    return OutputFormat::npz;
  else
    throw BlacklightException("Unknown string used for OutputFormat value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as Coordinates enumerations
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid Coordinates
// Notes:
//   Valid options:
//     "sph_ks": spherical Kerr-Schild
//     "cart_ks": Cartesian Kerr-Schild
Coordinates InputReader::ReadCoordinates(const std::string &string)
{
  if (string == "sph_ks")
    return Coordinates::sph_ks;
  else if (string == "cart_ks")
    return Coordinates::cart_ks;
  else
    throw BlacklightException("Unknown string used for Coordinates value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as Camera enumerations
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid Camera
// Notes:
//   Valid options:
//     "plane": rays originate from plane, parallel to central ray (which is perpendicular to
//         plane); appropriate for images as seen from infinity
//     "pinhole": rays originate from point, going in different directions; appropriate for camera
//         located near source
Camera InputReader::ReadCamera(const std::string &string)
{
  if (string == "plane")
    return Camera::plane;
  else if (string == "pinhole")
    return Camera::pinhole;
  else
    throw BlacklightException("Unknown string used for Camera value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as FrequencyNormalization enumerations
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid FrequencyNormalization
// Notes:
//   Valid options:
//     "camera": input image_frequency is taken to be the frequency as seen by the center of the
//         camera, accounting for its position and velocity; that is, the covariant time component
//         of photon momentum in the camera frame at the camera location is -image_frequency
//     "infinity": input image_frequency is taken to be the frequency of the light at the center of
//         the camera were it to be transported along geodesics to infinity and measured by an
//         observer at rest; that is, the covariant time component of photon momentum in the
//         coordinate frame is -image_frequency
FrequencyNormalization InputReader::ReadFrequencyNormalization(const std::string &string)
{
  if (string == "camera")
    return FrequencyNormalization::camera;
  else if (string == "infinity")
    return FrequencyNormalization::infinity;
  else
    throw BlacklightException("Unknown string used for FrequencyNormalization value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as RayTerminate enumerations
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid RayTerminate
// Notes:
//   Valid options:
//     "photon": rays are terminated upon reaching the prograde equatorial photon orbit radius
//     "multiplicative": rays are terminated upon reaching the horizon radius times ray_factor
//         (dimensionless)
//     "additive": rays are terminated upon reaching the horizon radius plus ray_factor
//         (gravitational units)
RayTerminate InputReader::ReadRayTerminate(const std::string &string)
{
  if (string == "photon")
    return RayTerminate::photon;
  else if (string == "multiplicative")
    return RayTerminate::multiplicative;
  else if (string == "additive")
    return RayTerminate::additive;
  else
    throw BlacklightException("Unknown string used for RayTerminate value.");
}
