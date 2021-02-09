// Blacklight output writer header

#ifndef WRITE_OUTPUT_H_
#define WRITE_OUTPUT_H_

// C++ headers
#include <cstdint>  // uint8_t, uint16_t, uint32_t
#include <cstdio>   // size_t
#include <fstream>  // ofstream
#include <string>   // string

// Blacklight headers
#include "array.hpp"       // Array
#include "blacklight.hpp"  // enumerations
#include "ray_tracer.hpp"  // RayTracer
#include "read_input.hpp"  // InputReader

//--------------------------------------------------------------------------------------------------

// Input reader
struct OutputWriter
{
  // Constructors and destructor
  OutputWriter(const InputReader *p_input_reader, const RayTracer *p_ray_tracer);
  OutputWriter(const OutputWriter &source) = delete;
  OutputWriter &operator=(const OutputWriter &source) = delete;
  ~OutputWriter() {}

  // Input data - output parameters
  OutputFormat output_format;
  std::string output_file;
  bool output_camera;

  // Input data - image parameters
  Camera im_camera;

  // File data
  std::ofstream *p_output_stream;

  // Image data
  Array<double> image;
  Array<double> im_pos;
  Array<double> im_dir;
  double im_width;
  double mass_msun;

  // External function
  void Write();

  // Internal functions
  void WriteRaw();
  void WriteNpy();
  void WriteNpz();
  std::size_t GenerateNpyFromDoubleArray(const Array<double> &array, uint8_t **p_buffer);
  std::size_t GenerateZIPLocalFileHeader(const uint8_t *record, std::size_t record_length,
      const char *record_name, uint8_t **p_buffer);
  std::size_t GenerateZIPCentralDirectoryHeader(const uint8_t *local_header,
      std::size_t local_header_offset, uint8_t **p_buffer);
  std::size_t GenerateZIPEndOfCentralDirectoryRecord(std::size_t central_directory_offset,
      std::size_t central_directory_length, int num_central_directory_entries, uint8_t **p_buffer);
  uint32_t CalculateCRC32(const uint8_t *message, std::size_t message_length);
};

#endif
