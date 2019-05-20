// Ray Trace input reader

// C++ headers
#include <algorithm>   // remove_if
#include <cctype>      // isspace
#include <fstream>     // ifstream
#include <iostream>    // cout
#include <string>      // getline, string

// Ray Trace headers
#include "read_input.hpp"
#include "exceptions.hpp"  // ray_trace_exception

// Input reader constructor
// Inputs:
//   input_file: C string containing name of input file

input_reader::input_reader(const std::string input_file_)
  : input_file(input_file_)
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
    if (line.empty() or line[1] == '#')
      continue;
  }
}

// Input reader destructor

input_reader::~input_reader() {}

// Definition of what constitutes a space
// Inputs:
//   c: character to be tested
// Outputs:
//   returned value: flag indicating character is a removable space

bool input_reader::removeable_space(unsigned char c)
{
  return std::isspace(c) != 0;
}
