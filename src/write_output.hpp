// Ray Trace output writer header

#ifndef WRITE_OUTPUT_H_
#define WRITE_OUTPUT_H_

// C++ headers
#include <string>  // string

// Ray Trace headers
#include "array.hpp"       // Array
#include "ray_tracer.hpp"  // RayTracer

//--------------------------------------------------------------------------------------------------

// Input reader
struct OutputWriter
{
  // Constructor
  OutputWriter(const std::string output_file_, const RayTracer &ray_tracer);

  // Metadata
  const std::string output_file;

  // Image data
  Array<double> image;

  // Functions
  void Write();
};

#endif
