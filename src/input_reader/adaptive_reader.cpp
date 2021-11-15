// Blacklight input reader - adaptive reader

// C++ headers
#include <optional>  // optional
#include <sstream>   // ostringstream
#include <string>    // stod, stoi, string

// Blacklight headers
#include "input_reader.hpp"
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function for parsing adaptive options
// Inputs:
//   key: input key as a string without leading "adaptive_" or "adaptive_region_"
//   val: input value as a string
// Outputs: (none)
// Notes:
//   Values are only recorded if space is allocated for them, so the keys should occur in a sensible
//       order in the input file.
//   Variables indexed beyond what is allocated (e.g. adaptive_region_2_level, when
//       adaptive_num_regions = 1) are silently ignored.
void InputReader::ReadAdaptive(const std::string &key, const std::string &val)
{
  // Read total number of forced refinement regions
  if (key == "num_regions")
  {
    adaptive_num_regions = std::stoi(val);
    adaptive_region_levels = new std::optional<int>[adaptive_num_regions.value()];
    adaptive_region_x_min_vals = new std::optional<double>[adaptive_num_regions.value()];
    adaptive_region_x_max_vals = new std::optional<double>[adaptive_num_regions.value()];
    adaptive_region_y_min_vals = new std::optional<double>[adaptive_num_regions.value()];
    adaptive_region_y_max_vals = new std::optional<double>[adaptive_num_regions.value()];
  }

  // Read minimum refinement level for a particular region
  else if (key.size() >= 7 and key.compare(key.size() - 6, key.npos, "_level") == 0)
  {
    int region_num = std::stoi(key.substr(0, key.size() - 6)) - 1;
    if (region_num >= adaptive_num_regions.value())
      return;
    adaptive_region_levels[region_num] = std::stoi(val);
  }

  // Read left boundary for a particular region
  else if (key.size() >= 7 and key.compare(key.size() - 6, key.npos, "_x_min") == 0)
  {
    int region_num = std::stoi(key.substr(0, key.size() - 6)) - 1;
    if (region_num >= adaptive_num_regions.value())
      return;
    adaptive_region_x_min_vals[region_num] = std::stod(val);
  }

  // Read right boundary for a particular region
  else if (key.size() >= 7 and key.compare(key.size() - 6, key.npos, "_x_max") == 0)
  {
    int region_num = std::stoi(key.substr(0, key.size() - 6)) - 1;
    if (region_num >= adaptive_num_regions.value())
      return;
    adaptive_region_x_max_vals[region_num] = std::stod(val);
  }

  // Read bottom boundary for a particular region
  else if (key.size() >= 7 and key.compare(key.size() - 6, key.npos, "_y_min") == 0)
  {
    int region_num = std::stoi(key.substr(0, key.size() - 6)) - 1;
    if (region_num >= adaptive_num_regions.value())
      return;
    adaptive_region_y_min_vals[region_num] = std::stod(val);
  }

  // Read top boundary for a particular region
  else if (key.size() >= 7 and key.compare(key.size() - 6, key.npos, "_y_max") == 0)
  {
    int region_num = std::stoi(key.substr(0, key.size() - 6)) - 1;
    if (region_num >= adaptive_num_regions.value())
      return;
    adaptive_region_y_max_vals[region_num] = std::stod(val);
  }

  // Handle unknown entry
  else
  {
    std::ostringstream message;
    message << "Unknown key (adaptive_region_" << key << ") in input file.";
    throw BlacklightException(message.str().c_str());
  }
  return;
}
