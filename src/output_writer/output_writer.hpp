// Blacklight output writer header

#ifndef OUTPUT_WRITER_H_
#define OUTPUT_WRITER_H_

// C++ headers
#include <cstddef>  // size_t
#include <cstdint>  // uint8_t, uint32_t
#include <fstream>  // ofstream
#include <string>   // string

// Blacklight headers
#include "../blacklight.hpp"                                 // enums
#include "../geodesic_integrator/geodesic_integrator.hpp"    // GeodesicIntegrator
#include "../input_reader/input_reader.hpp"                  // InputReader
#include "../radiation_integrator/radiation_integrator.hpp"  // RadiationIntegrator
#include "../utils/array.hpp"                                // Array

//--------------------------------------------------------------------------------------------------

// Output reader
struct OutputWriter
{
  // Constructors and destructor
  OutputWriter(const InputReader *p_input_reader_, const GeodesicIntegrator *p_geodesic_integrator_,
      const RadiationIntegrator *p_radiation_integrator_);
  OutputWriter(const OutputWriter &source) = delete;
  OutputWriter &operator=(const OutputWriter &source) = delete;
  ~OutputWriter();

  // Pointers to other objects
  const InputReader *p_input_reader;
  const GeodesicIntegrator *p_geodesic_integrator;
  const RadiationIntegrator *p_radiation_integrator;

  // Input data - general
  ModelType model_type;

  // Input data - output parameters
  OutputFormat output_format;
  std::string output_file;
  bool output_camera;

  // Input data - simulation parameters
  bool simulation_multiple;

  // Input data - camera parameters
  Camera camera_type;
  int camera_resolution;

  // Input data - image parameters
  bool image_light;
  bool image_polarization;
  bool image_time;
  bool image_length;
  bool image_lambda;
  bool image_emission;
  bool image_tau;
  bool image_lambda_ave;
  bool image_emission_ave;
  bool image_tau_int;
  int image_num_cell_values;
  int image_offset_time;
  int image_offset_length;
  int image_offset_lambda;
  int image_offset_emission;
  int image_offset_tau;
  int image_offset_lambda_ave;
  int image_offset_emission_ave;
  int image_offset_tau_int;

  // Input data - adaptive parameters
  bool adaptive_on;
  int adaptive_block_size;
  int adaptive_max_level;

  // Flag for tracking function calls
  bool first_time = true;

  // File data
  std::ofstream *p_output_stream;

  // Image data
  Array<double> image;
  Array<double> camera_pos;
  Array<double> camera_dir;
  Array<double> camera_width_array;
  Array<double> image_frequency_array;
  Array<double> mass_msun_array;

  // Adaptive data
  Array<int> adaptive_num_levels_array;
  Array<int> block_counts_array;
  Array<double> *image_adaptive;
  Array<int> *camera_loc_adaptive;
  Array<double> *camera_pos_adaptive;
  Array<double> *camera_dir_adaptive;

  // Name data
  static constexpr int num_cell_names = 7;
  const char *cell_names[num_cell_names] =
      {"rho", "n_e", "p_gas", "Theta_e", "B", "sigma", "beta_inverse"};

  // External function
  void Write(int snapshot);

  // Internal functions - output_writer.cpp
  std::string FormatFilename(int file_number);

  // Internal functions - raw_format.cpp
  void WriteRaw();

  // Internal functions - numpy_format.cpp
  void WriteNpy();
  void WriteNpz();
  std::size_t GenerateNpyFromArray(const Array<double> &array, int num_dims, uint8_t **p_buffer);
  std::size_t GenerateNpyFromArray(const Array<int> &array, int num_dims, uint8_t **p_buffer);

  // Internal functions - zip_format.cpp
  std::size_t GenerateZIPLocalFileHeader(const uint8_t *record, std::size_t record_length,
      const char *record_name, uint8_t **p_buffer);
  std::size_t GenerateZIPCentralDirectoryHeader(const uint8_t *local_header,
      std::size_t local_header_offset, uint8_t **p_buffer);
  std::size_t GenerateZIPEndOfCentralDirectoryRecord(std::size_t central_directory_offset,
      std::size_t central_directory_length, int num_central_directory_entries, uint8_t **p_buffer);
  uint32_t CalculateCRC32(const uint8_t *message, std::size_t message_length);
};

#endif
