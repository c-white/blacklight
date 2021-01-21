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
  // Constructors and destructor
  OutputWriter(const InputReader *p_input_reader, const RayTracer *p_ray_tracer);
  OutputWriter(const OutputWriter &source) = delete;
  OutputWriter &operator=(const OutputWriter &source) = delete;
  ~OutputWriter() {}

  // Input data - output parameters
  OutputFormat output_format;
  std::string output_file;

  // Image data
  Array<double> image;

  // Functions
  void Write();
};

#endif
