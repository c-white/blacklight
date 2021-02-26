// Blacklight output writer

// C++ headers
#include <fstream>   // ofstream
#include <ios>       // ios_base
#include <optional>  // optional

// Blacklight headers
#include "output_writer.hpp"
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../ray_tracer/ray_tracer.hpp"      // RayTracer
#include "../utils/array.hpp"                // Array
#include "../utils/exceptions.hpp"           // BlacklightException

//--------------------------------------------------------------------------------------------------

// Output writer constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
//   p_ray_tracer: pointer to object containing camera and processed image
// Notes:
//   File is not opened for writing until Write() function is called, in order to not bother keeping
//       an open file idle for the bulk of the computation.
OutputWriter::OutputWriter(const InputReader *p_input_reader, const RayTracer *p_ray_tracer)
{
  // Copy output parameters
  output_format = p_input_reader->output_format.value();
  output_file = p_input_reader->output_file_formatted;
  if (output_format == OutputFormat::npz)
  {
    output_params = p_input_reader->output_params.value();
    output_camera = p_input_reader->output_camera.value();
  }

  // Copy image parameters
  if (output_format == OutputFormat::npz and output_params)
  {
    image_width_array.Allocate(1);
    image_frequency_array.Allocate(1);
    mass_msun_array.Allocate(1);
    image_width_array(0) = p_input_reader->image_width.value();
    image_frequency_array(0) = p_input_reader->image_frequency.value();
    mass_msun_array(0) = p_ray_tracer->mass_msun;
  }

  // Make local shallow copies of data
  image = p_ray_tracer->image;
  if (output_format == OutputFormat::npz and output_camera)
    image_camera = p_input_reader->image_camera.value();
  if (output_format == OutputFormat::npz and output_camera and image_camera == Camera::plane)
    image_position = p_ray_tracer->image_position;
  if (output_format == OutputFormat::npz and output_camera and image_camera == Camera::pinhole)
    image_direction = p_ray_tracer->image_direction;
}

//--------------------------------------------------------------------------------------------------

// Output writer write function
// Inputs: (none)
// Outputs: (none)
void OutputWriter::Write()
{
  // Open output file
  p_output_stream = new std::ofstream(output_file, std::ios_base::out | std::ios_base::binary);
  if (not p_output_stream->is_open())
    throw BlacklightException("Could not open output file.");

  // Write image data based on desired file format
  switch (output_format)
  {
    case OutputFormat::raw:
    {
      WriteRaw();
      break;
    }
    case OutputFormat::npy:
    {
      WriteNpy();
      break;
    }
    case OutputFormat::npz:
    {
      WriteNpz();
      break;
    }
  }

  // Close output file
  delete p_output_stream;
  return;
}
