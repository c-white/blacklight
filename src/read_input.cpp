// Blacklight input reader

// C++ headers
#include <algorithm>   // remove_if
#include <cctype>      // isspace
#include <fstream>     // ifstream
#include <string>      // getline, stod, stoi, string

// Blacklight headers
#include "read_input.hpp"
#include "blacklight.hpp"  // math, Coordinates
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

    // Store file name data
    if (key == "data_file")
      data_file = val;
    else if (key == "output_file")
      output_file = val;

    // Store coordinate data
    else if (key == "bh_m")
      bh_m = std::stod(val);
    else if (key == "bh_a")
      bh_a = std::stod(val);
    else if (key == "coord")
      coord = ReadCoordinates(val);

    // Store unit data
    else if (key == "m_msun")
      m_msun = std::stod(val);
    else if (key == "rho_unit")
      rho_unit = std::stod(val);

    // Store plasma data
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

    // Store image data
    else if (key == "im_radius")
      im_r = std::stod(val);
    else if (key == "im_theta")
      im_th = ReadPole(val, &im_pole) * math::pi/180.0;
    else if (key == "im_phi")
      im_ph = std::stod(val) * math::pi/180.0;
    else if (key == "im_rot")
      im_rot = std::stod(val) * math::pi/180.0;
    else if (key == "im_width")
      im_width = std::stod(val);
    else if (key == "im_res")
      im_res = std::stoi(val);
    else if (key == "im_freq")
      im_freq = std::stod(val);

    // Store ray data
    else if (key == "ray_step")
      ray_step = std::stod(val);
    else if (key == "ray_max_steps")
      ray_max_steps = std::stoi(val);
    else if (key == "ray_flat")
      ray_flat = ReadBool(val);

    // Store performance data
    else if (key == "num_threads")
      num_threads = std::stoi(val);

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
//   p_pole_flag: set to true only if input string converts to 0.0 or 180.0
double InputReader::ReadPole(const std::string &string, bool *p_pole_flag)
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
    return sph_ks;
  else if (string == "cart_ks")
    return cart_ks;
  else
    throw BlacklightException("Unknown string used for Coordinates value.");
}
