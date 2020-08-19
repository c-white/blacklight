// Blacklight output writer header

#ifndef WRITE_OUTPUT_H_
#define WRITE_OUTPUT_H_

// C++ headers
#include <string>  // string

// Blacklight headers
#include "array.hpp"       // Array
#include "blacklight.hpp"  // enumerations
#include "ray_tracer.hpp"  // RayTracer
#include "read_input.hpp"  // InputReader

//--------------------------------------------------------------------------------------------------

// Input reader
struct OutputWriter
{
  // Constructor
  OutputWriter(const InputReader &input_reader, const RayTracer &ray_tracer);

  // Input data - output parameters
  OutputFormat output_format;
  std::string output_file;

  // Image data
  Array<double> image;

  // Functions
  void Write();
};

#endif
