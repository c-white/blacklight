// Ray Trace input reader

// C++ headers
#include <algorithm>   // remove_if
#include <cctype>      // isspace
#include <fstream>     // ifstream
#include <string>      // getline, stod, string

// Ray Trace headers
#include "read_input.hpp"
#include "exceptions.hpp"  // ray_trace_exception
#include "ray_trace.hpp"   // math

//--------------------------------------------------------------------------------------------------

// Input reader constructor
// Inputs:
//   input_file_: name of input file
input_reader::input_reader(const std::string input_file_)
  : input_file(input_file_) {}

//--------------------------------------------------------------------------------------------------

// Input reader read and initialize function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Initializes all member objects.
void input_reader::read()
{
  // Open input file
  std::ifstream input_stream(input_file);
  if (not input_stream.is_open())
    throw ray_trace_exception("Error: Could not open input file.\n");

  // Process file line by line
  for (std::string line; std::getline(input_stream, line); )
  {
    // Remove spaces
    line.erase(std::remove_if(line.begin(), line.end(), removeable_space), line.end());

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
      throw ray_trace_exception("Error: Invalid assignment in input file.\n");
    std::string key = line.substr(0, pos);
    std::string val = line.substr(pos + 1, line.size());

    // Store values
    if (key == "data_file")
      data_file = val;
    else if (key == "output_file")
      output_file = val;
    else if (key == "m")
      m = stod(val);
    else if (key == "a")
      a = stod(val);
    else if (key == "im_theta")
      im_th = stod(val) * math::pi/180.0;
    else if (key == "im_phi")
      im_ph = stod(val) * math::pi/180.0;
    else if (key == "im_width")
      im_width = stod(val);
    else if (key == "im_res")
      im_res = stoi(val);
    else if (key == "num_samples")
      num_samples = stoi(val);
    else
      throw ray_trace_exception("Error: Unknown key in input file.\n");
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Definition of what constitutes a space
// Inputs:
//   c: character to be tested
// Outputs:
//   returned value: flag indicating character is a removable space
bool input_reader::removeable_space(unsigned char c)
{
  return std::isspace(c) != 0;
}
