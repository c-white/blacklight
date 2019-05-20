// Ray Trace Athena++ reader

// C++ headers
#include <fstream>  // ifstream
#include <string>   // string

// Ray Trace headers
#include "read_athena.hpp"
#include "read_input.hpp"   // input_reader
#include "exceptions.hpp"   // ray_trace_exception

// Athena++ reader constructor
// Inputs:
//   input_file: object containing input file data

athena_reader::athena_reader(const std::string data_file_)
  : data_file(data_file_) {}

// Athena++ reader destructor

athena_reader::~athena_reader() {}

// Athena++ reader read and initialize function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Initializes all member objects

void athena_reader::read()
{
  // Open data file
  std::ifstream data_stream(data_file);
  if (not data_stream.is_open())
    throw ray_trace_exception("Error: Could not open data file.\n");
}
