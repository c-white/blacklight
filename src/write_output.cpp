// Blacklight output writer

// C++ headers
#include <fstream>  // ofstream
#include <ios>      // ios_base, streamsize
#include <string>   // string

// Blacklight headers
#include "write_output.hpp"
#include "array.hpp"         // Array
#include "blacklight.hpp"    // enumerations
#include "exceptions.hpp"    // BlacklightException
#include "ray_tracer.hpp"    // RayTracer
#include "read_input.hpp"    // InputReader

//--------------------------------------------------------------------------------------------------

// Output writer constructor
// Inputs:
//   output_file_: name of output file
//   ray_tracer: object containing processed image
OutputWriter::OutputWriter(const InputReader &input_reader, const RayTracer &ray_tracer)
{
  // Copy output parameters
  output_format = input_reader.output_format;
  output_file = input_reader.output_file;

  // Make local shallow copy of image data
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

  // Write image data based on desired file format
  switch (output_format)
  {
    // Raw data without metadata
    case raw:
    default:
    {
      const char *data_pointer = reinterpret_cast<char *>(image.data);
      std::streamsize data_size = static_cast<std::streamsize>(image.GetNumBytes());
      output_stream.write(data_pointer, data_size);
      break;
    }
  }
  return;
}
