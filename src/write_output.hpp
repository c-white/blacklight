// Ray Trace output writer header

#ifndef WRITE_OUTPUT_H_
#define WRITE_OUTPUT_H_

// C++ headers
#include <ios>     // streamsize
#include <string>  // string

// Ray Trace headers
#include "ray_tracer.hpp"  // RayTracer

//--------------------------------------------------------------------------------------------------

// Input reader
struct OutputWriter
{
  // Constructor
  OutputWriter(const std::string output_file_, const RayTracer &ray_tracer);

  // Metadata
  const std::string output_file;
  const char * const data_pointer;
  const std::streamsize data_size;

  // Functions
  void Write();
};

#endif
