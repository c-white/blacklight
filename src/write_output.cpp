// Blacklight output writer

// C++ headers
#include <cstdio>   // size_t, snprintf
#include <cstring>  // memcpy, memset
#include <fstream>  // ofstream
#include <ios>      // ios_base, streamsize

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
// Notes:
//   Must be run on little-endian machine.
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
    {
      const char *data_pointer = reinterpret_cast<char *>(image.data);
      std::streamsize data_size = static_cast<std::streamsize>(image.GetNumBytes());
      output_stream.write(data_pointer, data_size);
      break;
    }

    // Numpy file with single array and minimal metadata
    case npy:
    {
      // Prepare buffer for header
      const std::size_t buffer_length = 128;
      char header_buffer[buffer_length];
      std::memset(header_buffer, ' ', buffer_length - 1);
      header_buffer[buffer_length-1] = '\n';

      // Write magic string and version number to buffer
      std::size_t length = 0;
      const char *magic_string = "\x93NUMPY";
      std::memcpy(header_buffer, magic_string, 6);
      length += 6;
      const char *version = "\x01\x00";
      std::memcpy(header_buffer + length, version, 2);
      length += 2;

      // Write length of header proper to buffer
      unsigned short int header_len = static_cast<unsigned short int>(buffer_length - length - 2);
      std::memcpy(header_buffer + length, reinterpret_cast<const char *>(&header_len), 2);
      length += 2;

      // Write header proper to buffer
      int num_written = std::snprintf(header_buffer + length, buffer_length - length - 1,
          "{'descr': '<f8', 'fortran_order': False, 'shape': (%d,%d)}", image.n2, image.n1);
      if (num_written < 0)
        throw BlacklightException("Error writing output file.");
      length += static_cast<std::size_t>(num_written);
      header_buffer[length] = ' ';

      // Write buffer
      output_stream.write(header_buffer, buffer_length);

      // Write data
      const char *data_pointer = reinterpret_cast<char *>(image.data);
      std::streamsize data_size = static_cast<std::streamsize>(image.GetNumBytes());
      output_stream.write(data_pointer, data_size);
      break;
    }
  }
  return;
}
