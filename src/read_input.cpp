// Ray Trace input reader

// C++ headers
#include <algorithm>   // remove_if
#include <cctype>      // isspace
#include <fstream>     // ifstream
#include <string>      // getline, stod, string

// Ray Trace headers
#include "read_input.hpp"
#include "exceptions.hpp"  // ray_trace_exception

//--------------------------------------------------------------------------------------------------

// Input reader constructor
// Inputs:
//   input_file_: name of input file
input_reader::input_reader(const std::string input_file_)
  : input_file(input_file_) {}

//--------------------------------------------------------------------------------------------------

// Input reader destructor
input_reader::~input_reader() {}

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

    // Skip blank lines and comments
    if (line.empty() or line[0] == '#')
      continue;

    // Split on '='
    std::string::size_type pos = line.find_first_of("=");
    if (pos == std::string::npos)
      throw ray_trace_exception("Error: Invalid assignment in input file.\n");
    std::string key = line.substr(0, pos);
    std::string val = line.substr(pos + 1, line.size());

    // Store values
    if (key == "data_file")
      data_file = val;
    else if (key == "m")
      m = stod(val);
    else if (key == "a")
      a = stod(val);
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
