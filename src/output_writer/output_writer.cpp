// Blacklight output writer

// C++ headers
#include <fstream>   // ofstream
#include <ios>       // ios_base
#include <optional>  // optional

// Blacklight headers
#include "output_writer.hpp"
#include "../blacklight.hpp"                                 // enums
#include "../geodesic_integrator/geodesic_integrator.hpp"    // GeodesicIntegrator
#include "../input_reader/input_reader.hpp"                  // InputReader
#include "../radiation_integrator/radiation_integrator.hpp"  // RadiationIntegrator
#include "../utils/array.hpp"                                // Array
#include "../utils/exceptions.hpp"                           // BlacklightException

//--------------------------------------------------------------------------------------------------

// Output writer constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
//   p_geodesic_integrator: pointer to object containing ray data
//   p_radiation_integrator: pointer to object containing processed image
// Notes:
//   File is not opened for writing until Write() function is called, in order to not bother keeping
//       an open file idle for the bulk of the computation, and because the file name might be
//       reformatted after this constructor is called.
OutputWriter::OutputWriter(const InputReader *p_input_reader_,
    const GeodesicIntegrator *p_geodesic_integrator,
    const RadiationIntegrator *p_radiation_integrator_)
  : p_input_reader(p_input_reader_), p_radiation_integrator(p_radiation_integrator_)
{
  // Copy output parameters
  output_format = p_input_reader->output_format.value();
  if (output_format == OutputFormat::npz)
  {
    output_params = p_input_reader->output_params.value();
    output_camera = p_input_reader->output_camera.value();
  }

  // Copy image parameters
  image_resolution = p_input_reader->image_resolution.value();
  if (output_format == OutputFormat::npz and output_params)
  {
    image_width_array.Allocate(1);
    image_frequency_array.Allocate(1);
    mass_msun_array.Allocate(1);
    image_width_array(0) = p_input_reader->image_width.value();
    image_frequency_array(0) = p_input_reader->image_frequency.value();
    mass_msun_array(0) = p_radiation_integrator->mass_msun;
  }

  // Make local shallow copies of data, reshaping the arrays
  if (output_format == OutputFormat::npz and output_camera)
    image_camera = p_input_reader->image_camera.value();
  if (output_format == OutputFormat::npz and output_camera and image_camera == Camera::plane)
  {
    camera_pos = p_geodesic_integrator->camera_pos;
    camera_pos.n3 = image_resolution;
    camera_pos.n2 = image_resolution;
  }
  if (output_format == OutputFormat::npz and output_camera and image_camera == Camera::pinhole)
  {
    camera_dir = p_geodesic_integrator->camera_dir;
    camera_dir.n3 = image_resolution;
    camera_dir.n2 = image_resolution;
  }
}

//--------------------------------------------------------------------------------------------------

// Output writer write function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Opens and closes stream for reading.
void OutputWriter::Write()
{
  // Make shallow copy of image
  if (first_time)
  {
    image = p_radiation_integrator->image;
    image.n3 = image.n2;
    image.n2 = image_resolution;
    image.n1 = image_resolution;
  }

  // Open output file
  output_file = p_input_reader->output_file_formatted;
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

  // Update first time flag
  first_time = false;
  return;
}
