// Ray Trace Athena++ reader

// C++ headers
#include <fstream>  // ifstream

// Ray Trace headers
#include "read_athena.hpp"
#include "read_input.hpp"   // input_reader
#include "exceptions.hpp"   // ray_trace_exception

// Athena++ reader constructor
// Inputs:
//   input_file: object containing input file data

athena_reader::athena_reader(const input_reader &input_reader)
{
  // Open data file
  std::ifstream data_stream(input_reader.data_file);
  if (not data_stream.is_open())
    throw ray_trace_exception("Error: Could not open data file.\n");
}

// Athena++ reader destructor

athena_reader::~athena_reader() {}
