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
    simulation_multiple = p_input_reader->simulation_multiple.value();

  // Copy camera parameters
  if (output_format == OutputFormat::npz and output_camera)
    camera_type = p_input_reader->camera_type.value();
  camera_resolution = p_input_reader->camera_resolution.value();

  // Copy image parameters
  image_light = p_input_reader->image_light.value();
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
  image_num_cell_values = p_radiation_integrator->image_num_cell_values;
  if (image_num_cell_values != num_cell_names)
    throw BlacklightException("Unexpected number of cell values.");
  image_offset_time = p_radiation_integrator->image_offset_time;
  image_offset_length = p_radiation_integrator->image_offset_length;
  image_offset_lambda = p_radiation_integrator->image_offset_lambda;
  image_offset_emission = p_radiation_integrator->image_offset_emission;
  image_offset_tau = p_radiation_integrator->image_offset_tau;
  image_offset_lambda_ave = p_radiation_integrator->image_offset_lambda_ave;
  image_offset_emission_ave = p_radiation_integrator->image_offset_emission_ave;
  image_offset_tau_int = p_radiation_integrator->image_offset_tau_int;

  // Copy adaptive parameters
  adaptive_on = p_input_reader->adaptive_on.value();
  if (adaptive_on)
  {
    if (output_format != OutputFormat::npz)
      throw BlacklightException("Only npz outputs support adaptive ray tracing.");
    adaptive_block_size = p_input_reader->adaptive_block_size.value();
    adaptive_max_level = p_input_reader->adaptive_max_level.value();
  }

  // Copy scalar parameters to arrays
  if (output_format == OutputFormat::npz)
  {
    camera_width_array.Allocate(1);
    image_frequency_array.Allocate(1);
    mass_msun_array.Allocate(1);
    camera_width_array(0) = p_input_reader->camera_width.value();
    image_frequency_array(0) = p_input_reader->image_frequency.value();
    mass_msun_array(0) = p_radiation_integrator->mass_msun;
  }

  // Make shallow copies of camera data, reshaping the arrays
  if (output_format == OutputFormat::npz and output_camera and camera_type == Camera::plane)
  {
    camera_pos = p_geodesic_integrator->camera_pos;
    camera_pos.n3 = camera_resolution;
    camera_pos.n2 = camera_resolution;
  }
  if (output_format == OutputFormat::npz and output_camera and camera_type == Camera::pinhole)
  {
    camera_dir = p_geodesic_integrator->camera_dir;
    camera_dir.n3 = camera_resolution;
    camera_dir.n2 = camera_resolution;
  }

  // Allocate space for shallow copies of adaptive data
  adaptive_num_levels_array.Allocate(1);
  if (adaptive_on)
  {
    image_adaptive = new Array<double>[adaptive_max_level+1];
    camera_loc_adaptive = new Array<int>[adaptive_max_level+1];
    camera_pos_adaptive = new Array<double>[adaptive_max_level+1];
    camera_dir_adaptive = new Array<double>[adaptive_max_level+1];
  }
}

//--------------------------------------------------------------------------------------------------

// Output writer destructor
OutputWriter::~OutputWriter()
{
  if (adaptive_on)
  {
    for (int level = 0; level <= adaptive_max_level; level++)
    {
      image_adaptive[level].Deallocate();
      camera_loc_adaptive[level].Deallocate();
      camera_pos_adaptive[level].Deallocate();
      camera_dir_adaptive[level].Deallocate();
    }
    delete[] image_adaptive;
    delete[] camera_loc_adaptive;
    delete[] camera_pos_adaptive;
    delete[] camera_dir_adaptive;
  }
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
  // Make shallow copy of root image data
  if (first_time)
  {
    image = p_radiation_integrator->image;
    image.n3 = image.n2;
    image.n2 = camera_resolution;
    image.n1 = camera_resolution;
  }

  // Make shallow copy of adaptive image data
  adaptive_num_levels_array(0) = p_radiation_integrator->adaptive_num_levels;
  if (adaptive_on)
  {
    block_counts_array.Allocate(adaptive_num_levels_array(0) + 1);
    block_counts_array(0) = p_radiation_integrator->block_counts[0];
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      block_counts_array(level) = p_radiation_integrator->block_counts[level];
      image_adaptive[level] = p_radiation_integrator->image_adaptive[level];
      image_adaptive[level].n4 = image_adaptive[level].n2;
      image_adaptive[level].n3 = block_counts_array(level);
      image_adaptive[level].n2 = adaptive_block_size;
      image_adaptive[level].n1 = adaptive_block_size;
    }
  }

  // Make shallow copy of adaptive camera data
  if (adaptive_on)
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
      camera_loc_adaptive[level] = p_geodesic_integrator->camera_loc_adaptive[level];
  if (adaptive_on and output_camera and camera_type == Camera::plane)
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      camera_pos_adaptive[level] = p_geodesic_integrator->camera_pos_adaptive[level];
      camera_pos_adaptive[level].n4 = block_counts_array(level);
      camera_pos_adaptive[level].n3 = adaptive_block_size;
      camera_pos_adaptive[level].n2 = adaptive_block_size;
    }
  if (adaptive_on and output_camera and camera_type == Camera::pinhole)
    for (int level = 1; level <= adaptive_num_levels_array(0); level++)
    {
      camera_dir_adaptive[level] = p_geodesic_integrator->camera_dir_adaptive[level];
      camera_dir_adaptive[level].n4 = block_counts_array(level);
      camera_dir_adaptive[level].n3 = adaptive_block_size;
      camera_dir_adaptive[level].n2 = adaptive_block_size;
    }

  // Open output file
  std::string output_file_formatted = output_file;
  if (model_type == ModelType::simulation and simulation_multiple)
    output_file_formatted = FormatFilename(snapshot);
  p_output_stream =
      new std::ofstream(output_file_formatted, std::ios_base::out | std::ios_base::binary);
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
