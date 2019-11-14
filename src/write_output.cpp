// Ray Trace output writer

// C++ headers
#include <fstream>  // ofstream
#include <ios>      // ios_base, streamsize
#include <string>   // string

// Ray Trace headers
#include "write_output.hpp"
#include "array.hpp"         // Array
#include "exceptions.hpp"    // RayTraceException
#include "ray_tracer.hpp"    // RayTracer

//--------------------------------------------------------------------------------------------------

// Output writer constructor
// Inputs:
//   output_file_: name of output file
//   ray_tracer: object containing processed image
OutputWriter::OutputWriter(const std::string output_file_, const RayTracer &ray_tracer)
  : output_file(output_file_),
    data_pointer(reinterpret_cast<char *>(ray_tracer.image.data)),
    data_size(static_cast<std::streamsize>(ray_tracer.image.GetNumBytes())) {}

//--------------------------------------------------------------------------------------------------

// Output writer write function
// Inputs: (none)
// Outputs: (none)
void OutputWriter::Write()
{
  // Open output file
  std::ofstream output_stream(output_file, std::ios_base::out | std::ios_base::binary);
  if (not output_stream.is_open())
    throw RayTraceException("Could not open output file.");

  // Write image data
  output_stream.write(data_pointer, data_size);
  return;
}
