// Blacklight output writer

// C++ headers
#include <fstream>  // ofstream
#include <ios>      // ios_base, streamsize
#include <string>   // string

// Blacklight headers
#include "write_output.hpp"
#include "array.hpp"         // Array
#include "exceptions.hpp"    // BlacklightException
#include "ray_tracer.hpp"    // RayTracer

//--------------------------------------------------------------------------------------------------

// Output writer constructor
// Inputs:
//   output_file_: name of output file
//   ray_tracer: object containing processed image
OutputWriter::OutputWriter(const std::string output_file_, const RayTracer &ray_tracer)
  : output_file(output_file_)
{
  image = ray_tracer.image;
}

//--------------------------------------------------------------------------------------------------

// Output writer write function
// Inputs: (none)
// Outputs: (none)
void OutputWriter::Write()
{
  // Open output file
  std::ofstream output_stream(output_file, std::ios_base::out | std::ios_base::binary);
  if (not output_stream.is_open())
    throw BlacklightException("Could not open output file.");

  // Write image data
  const char *data_pointer = reinterpret_cast<char *>(image.data);
  std::streamsize data_size = static_cast<std::streamsize>(image.GetNumBytes());
  output_stream.write(data_pointer, data_size);
  return;
}
