// Blacklight output writer

// C++ headers
#include <cstdio>    // snprintf
#include <fstream>   // ofstream
#include <ios>       // ios_base
#include <optional>  // optional
#include <sstream>   // ostringstream
#include <string>    // stoi, string

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
//   p_input_reader_: pointer to object containing input parameters
//   p_geodesic_integrator_: pointer to object containing ray data
//   p_radiation_integrator_: pointer to object containing processed image
// Notes:
//   File is not opened for writing until Write() function is called, in order to not bother keeping
//       an open file idle for the bulk of the computation, and because the file name might be
//       reformatted after this constructor is called.
OutputWriter::OutputWriter(const InputReader *p_input_reader_,
    const GeodesicIntegrator *p_geodesic_integrator_,
    const RadiationIntegrator *p_radiation_integrator_)
  : p_input_reader(p_input_reader_), p_geodesic_integrator(p_geodesic_integrator_),
    p_radiation_integrator(p_radiation_integrator_)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();

  // Copy output parameters
  output_format = p_input_reader->output_format.value();
  output_file = p_input_reader->output_file.value();
  if (output_format == OutputFormat::npz)
    output_camera = p_input_reader->output_camera.value();

  // Copy simulation parameters
  if (model_type == ModelType::simulation)
  {
    simulation_multiple = p_input_reader->simulation_multiple.value();
    if (simulation_multiple)
      simulation_start = p_input_reader->simulation_start.value();
  }

  // Copy camera parameters
  if (output_format == OutputFormat::npz and output_camera)
    camera_type = p_input_reader->camera_type.value();
  camera_resolution = p_input_reader->camera_resolution.value();

  // Copy image parameters
  image_light = p_input_reader->image_light.value();
  image_num_frequencies = p_input_reader->image_num_frequencies.value();
  if (image_num_frequencies > 1 and not (output_format == OutputFormat::npz))
    throw BlacklightException("Only npz support multiple frequencies.");
  if (image_light and model_type == ModelType::simulation)
  {
    image_polarization = p_input_reader->image_polarization.value();
    if (image_polarization
        and not (output_format == OutputFormat::npz or output_format == OutputFormat::npy))
      throw BlacklightException("Only npz or npy outputs support polarization.");
  }
  image_time = p_input_reader->image_time.value();
  image_length = p_input_reader->image_length.value();
  image_lambda = p_input_reader->image_lambda.value();
  image_emission = p_input_reader->image_emission.value();
  image_tau = p_input_reader->image_tau.value();
  image_lambda_ave = p_input_reader->image_lambda_ave.value();
  image_emission_ave = p_input_reader->image_emission_ave.value();
  image_tau_int = p_input_reader->image_tau_int.value();
  if ((image_time or image_length or image_lambda or image_emission or image_tau or image_lambda_ave
      or image_emission_ave or image_tau_int) and not (output_format == OutputFormat::npz))
    throw BlacklightException("Only npz outputs support non-light images.");
  image_offset_time = p_radiation_integrator->image_offset_time;
  image_offset_length = p_radiation_integrator->image_offset_length;
  image_offset_lambda = p_radiation_integrator->image_offset_lambda;
  image_offset_emission = p_radiation_integrator->image_offset_emission;
  image_offset_tau = p_radiation_integrator->image_offset_tau;
  image_offset_lambda_ave = p_radiation_integrator->image_offset_lambda_ave;
  image_offset_emission_ave = p_radiation_integrator->image_offset_emission_ave;
  image_offset_tau_int = p_radiation_integrator->image_offset_tau_int;

  // Copy rendering parameters
  render_num_images = p_input_reader->render_num_images.value();
  if (render_num_images > 0 and output_format != OutputFormat::npz)
    throw BlacklightException("Only npz outputs support rendering.");

  // Copy slow-light parameters
  if (model_type == ModelType::simulation)
  {
    slow_light_on = p_input_reader->slow_light_on.value();
    if (slow_light_on)
      slow_offset = p_input_reader->slow_offset.value();
  }

  // Copy adaptive parameters
  adaptive_max_level = p_input_reader->adaptive_max_level.value();
  if (adaptive_max_level > 0)
  {
    if (output_format != OutputFormat::npz)
      throw BlacklightException("Only npz outputs support adaptive ray tracing.");
    adaptive_block_size = p_input_reader->adaptive_block_size.value();
  }

  // Copy metadata arrays
  if (output_format == OutputFormat::npz)
  {
    mass_msun_array.Allocate(1);
    mass_msun_array(0) = p_radiation_integrator->mass_msun;
    camera_width_array.Allocate(1);
    camera_width_array(0) = p_input_reader->camera_width.value();
    image_frequencies = p_geodesic_integrator->image_frequencies;
  }

  // Allocate space for camera data
  camera_loc = new Array<int>[adaptive_max_level+1];
  camera_pos = new Array<double>[adaptive_max_level+1];
  camera_dir = new Array<double>[adaptive_max_level+1];

  // Allocate space for image data
  image = new Array<double>[adaptive_max_level+1];

  // Allocate space for render data
  render = new Array<double>[adaptive_max_level+1];

  // Allocate space for adaptive data
  adaptive_num_levels_array.Allocate(1);
}

//--------------------------------------------------------------------------------------------------

// Output writer destructor
OutputWriter::~OutputWriter()
{
  for (int level = 0; level <= adaptive_max_level; level++)
  {
    camera_loc[level].Deallocate();
    camera_pos[level].Deallocate();
    camera_dir[level].Deallocate();
    image[level].Deallocate();
    render[level].Deallocate();
  }
  delete[] camera_loc;
  delete[] camera_pos;
  delete[] camera_dir;
  delete[] image;
  delete[] render;
}

//--------------------------------------------------------------------------------------------------

// Output writer write function
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be written
// Outputs: (none)
// Notes:
//   Opens and closes stream for reading.
void OutputWriter::Write(int snapshot)
{
  // Copy adaptive data
  adaptive_num_levels_array(0) = p_radiation_integrator->adaptive_num_levels;
  if (adaptive_max_level > 0)
  {
    block_counts_array.Allocate(adaptive_num_levels_array(0) + 1);
    block_counts_array(0) = p_radiation_integrator->block_counts[0];
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
      block_counts_array(level) = p_radiation_integrator->block_counts[level];
  }

  // Make shallow copies of camera data, reshaping the arrays
  for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    camera_loc[level] = p_geodesic_integrator->camera_loc[level];
  if (output_format == OutputFormat::npz and output_camera and camera_type == Camera::plane)
  {
    camera_pos[0] = p_geodesic_integrator->camera_pos[0];
    camera_pos[0].n3 = camera_resolution;
    camera_pos[0].n2 = camera_resolution;
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      camera_pos[level] = p_geodesic_integrator->camera_pos[level];
      camera_pos[level].n4 = block_counts_array(level);
      camera_pos[level].n3 = adaptive_block_size;
      camera_pos[level].n2 = adaptive_block_size;
    }
  }
  if (output_format == OutputFormat::npz and output_camera and camera_type == Camera::pinhole)
  {
    camera_dir[0] = p_geodesic_integrator->camera_dir[0];
    camera_dir[0].n3 = camera_resolution;
    camera_dir[0].n2 = camera_resolution;
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      camera_dir[level] = p_geodesic_integrator->camera_dir[level];
      camera_dir[level].n4 = block_counts_array(level);
      camera_dir[level].n3 = adaptive_block_size;
      camera_dir[level].n2 = adaptive_block_size;
    }
  }

  // Make shallow copies of image data, reshaping the arrays
  if (image_light or image_time or image_length or image_lambda or image_emission or image_tau
      or image_lambda_ave or image_emission_ave or image_tau_int)
  {
    image[0] = p_radiation_integrator->image[0];
    image[0].n3 = image[0].n2;
    image[0].n2 = camera_resolution;
    image[0].n1 = camera_resolution;
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      image[level] = p_radiation_integrator->image[level];
      image[level].n4 = image[level].n2;
      image[level].n3 = block_counts_array(level);
      image[level].n2 = adaptive_block_size;
      image[level].n1 = adaptive_block_size;
    }
  }

  // Make shallow copies of render data, reshaping the arrays
  if (render_num_images > 0)
  {
    render[0] = p_radiation_integrator->render[0];
    render[0].n4 = render[0].n3;
    render[0].n3 = render[0].n2;
    render[0].n2 = camera_resolution;
    render[0].n1 = camera_resolution;
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      render[level] = p_radiation_integrator->render[level];
      render[level].n5 = render[level].n3;
      render[level].n4 = render[level].n2;
      render[level].n3 = block_counts_array(level);
      render[level].n2 = adaptive_block_size;
      render[level].n1 = adaptive_block_size;
    }
  }

  // Open output file
  std::string output_file_formatted = output_file;
  if (model_type == ModelType::simulation and simulation_multiple)
  {
    int file_number = snapshot + (slow_light_on ? slow_offset : simulation_start);
    output_file_formatted = FormatFilename(file_number);
  }
  p_output_stream =
      new std::ofstream(output_file_formatted, std::ios_base::out | std::ios_base::binary);
  if (not p_output_stream->is_open())
    throw BlacklightException("Could not open output file.");

  // Write image data based on desired file format
  if (output_format == OutputFormat::npz)
    WriteNpz();
  else if (output_format == OutputFormat::npy)
    WriteNpy();
  else if (output_format == OutputFormat::raw)
    WriteRaw();

  // Close output file
  delete p_output_stream;

  // Free memory
  block_counts_array.Deallocate();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to construct filename formatted with file number
// Inputs:
//   file_number: number of output file to construct
// Outputs:
//   returned_value: string containing formatted filename
std::string OutputWriter::FormatFilename(int file_number)
{
  // Locate braces
  std::string::size_type output_pos_open = output_file.find_first_of('{');
  if (output_pos_open == std::string::npos)
    throw BlacklightException("Invalid output_file for multiple runs.");
  std::string::size_type output_pos_close = output_file.find_first_of('}', output_pos_open);
  if (output_pos_close == std::string::npos)
    throw BlacklightException("Invalid output_file for multiple runs.");

  // Parse integer format string
  if (output_file[output_pos_close-1] != 'd')
    throw BlacklightException("Invalid output_file for multiple runs.");
  int output_field_length = 0;
  if (output_pos_close - output_pos_open > 2)
    output_field_length = std::stoi(output_file.substr(output_pos_open + 1,
        output_pos_close - output_pos_open - 2));
  int file_number_length = std::snprintf(nullptr, 0, "%d", file_number);
  if (file_number_length < 0)
    throw BlacklightException("Could not format file name.");
  int num_zeros = 0;
  if (file_number_length < output_field_length)
    num_zeros = output_field_length - file_number_length;

  // Create filename
  std::ostringstream output_filename;
  output_filename << output_file.substr(0, output_pos_open);
  for (int n = 0; n < num_zeros; n++)
    output_filename << "0";
  output_filename << file_number;
  output_filename << output_file.substr(output_pos_close + 1);
  std::string output_file_formatted = output_filename.str();
  return output_file_formatted;
}
