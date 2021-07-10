// Blacklight output writer - NumPy output formats

// C++ headers
#include <cstdint>  // uint8_t, uint16_t
#include <cstdio>   // size_t, snprintf
#include <cstring>  // memcpy, memset
#include <fstream>  // ofstream
#include <ios>      // streamsize

// Blacklight headers
#include "output_writer.hpp"
#include "../blacklight.hpp"        // enums
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Output writer for a NumPy .npy file
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
void OutputWriter::WriteNpy()
{
  uint8_t *image_buffer;
  std::size_t image_length = GenerateNpyFromDoubleArray(image, &image_buffer);
  char *buffer = reinterpret_cast<char *>(image_buffer);
  std::streamsize buffer_length = static_cast<std::streamsize>(image_length);
  p_output_stream->write(buffer, buffer_length);
  delete[] image_buffer;
  return;
}

//--------------------------------------------------------------------------------------------------

// Output writer for a NumPy .npz file
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
//   NumPy .npz files are simply ZIP files with each record containing a .npy file and named
//       according to the name of the array.
//   Implements ZIP format version 2.0: pkware.cachefly.net/webdocs/casestudies/APPNOTE.TXT
//   Does not implement ZIP64, and thus has a reasonably low limit to data sizes that should still
//       be plenty for almost all ray tracing applications.
void OutputWriter::WriteNpz()
{
  // Prepare buffers for data and headers
  int num_arrays = 1 + (output_params ? 3 : 0) + (output_camera ? 1 : 0);
  int camera_offset = 1 + (output_params ? 3 : 0);
  uint8_t **data_buffers = new uint8_t *[num_arrays];
  uint8_t **local_header_buffers = new uint8_t *[num_arrays];
  uint8_t **central_header_buffers = new uint8_t *[num_arrays];
  std::size_t *data_lengths = new std::size_t[num_arrays];
  std::size_t *local_header_lengths = new std::size_t[num_arrays];
  std::size_t *central_header_lengths = new std::size_t[num_arrays];

  // Write data as .npy files to buffers
  data_lengths[0] = GenerateNpyFromDoubleArray(image, &data_buffers[0]);
  if (output_params)
  {
    data_lengths[1] = GenerateNpyFromDoubleArray(image_width_array, &data_buffers[1]);
    data_lengths[2] = GenerateNpyFromDoubleArray(image_frequency_array, &data_buffers[2]);
    data_lengths[3] = GenerateNpyFromDoubleArray(mass_msun_array, &data_buffers[3]);
  }
  if (output_camera)
  {
    if (image_camera == Camera::plane)
      data_lengths[camera_offset] =
          GenerateNpyFromDoubleArray(image_position, &data_buffers[camera_offset]);
    else if (image_camera == Camera::pinhole)
      data_lengths[camera_offset] =
          GenerateNpyFromDoubleArray(image_direction, &data_buffers[camera_offset]);
  }

  // Write local file headers to buffers
  local_header_lengths[0] = GenerateZIPLocalFileHeader(data_buffers[0], data_lengths[0], "image",
      &local_header_buffers[0]);
  if (output_params)
  {
    local_header_lengths[1] = GenerateZIPLocalFileHeader(data_buffers[1], data_lengths[1],
        "image_width", &local_header_buffers[1]);
    local_header_lengths[2] = GenerateZIPLocalFileHeader(data_buffers[2], data_lengths[2],
        "image_frequency", &local_header_buffers[2]);
    local_header_lengths[3] = GenerateZIPLocalFileHeader(data_buffers[3], data_lengths[3],
        "mass_msun", &local_header_buffers[3]);
  }
  if (output_camera)
  {
    if (image_camera == Camera::plane)
      local_header_lengths[camera_offset] = GenerateZIPLocalFileHeader(data_buffers[camera_offset],
          data_lengths[camera_offset], "image_position", &local_header_buffers[camera_offset]);
    else if (image_camera == Camera::pinhole)
      local_header_lengths[camera_offset] = GenerateZIPLocalFileHeader(data_buffers[camera_offset],
          data_lengths[camera_offset], "image_direction", &local_header_buffers[camera_offset]);
  }

  // Write central directory headers to buffers
  std::size_t offset = 0;
  for (int n = 0; n < num_arrays; n++)
  {
    central_header_lengths[n] = GenerateZIPCentralDirectoryHeader(local_header_buffers[n], offset,
        &central_header_buffers[n]);
    offset += local_header_lengths[n] + data_lengths[n];
  }

  // Write end of central directory record to buffer
  std::size_t central_directory_length = 0;
  for (int n = 0; n < num_arrays; n++)
    central_directory_length += central_header_lengths[n];
  uint8_t *end_of_directory_buffer = nullptr;
  std::size_t end_of_directory_length = GenerateZIPEndOfCentralDirectoryRecord(offset,
      central_directory_length, num_arrays, &end_of_directory_buffer);

  // Write buffers to file
  for (int n = 0; n < num_arrays; n++)
  {
    char *buffer = reinterpret_cast<char *>(local_header_buffers[n]);
    std::streamsize length = static_cast<std::streamsize>(local_header_lengths[n]);
    p_output_stream->write(buffer, length);
    buffer = reinterpret_cast<char *>(data_buffers[n]);
    length = static_cast<std::streamsize>(data_lengths[n]);
    p_output_stream->write(buffer, length);
  }
  for (int n = 0; n < num_arrays; n++)
  {
    char *buffer = reinterpret_cast<char *>(central_header_buffers[n]);
    std::streamsize length = static_cast<std::streamsize>(central_header_lengths[n]);
    p_output_stream->write(buffer, length);
  }
  char *buffer = reinterpret_cast<char *>(end_of_directory_buffer);
  std::streamsize length = static_cast<std::streamsize>(end_of_directory_length);
  p_output_stream->write(buffer, length);

  // Free memory
  for (int n = 0; n < num_arrays; n++)
  {
    delete[] data_buffers[n];
    delete[] local_header_buffers[n];
    delete[] central_header_buffers[n];
  }
  delete[] data_buffers;
  delete[] local_header_buffers;
  delete[] central_header_buffers;
  delete[] end_of_directory_buffer;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for populating a buffer with the contents of a NumPy .npy file from a double Array
// Inputs:
//   array: Array to be written to buffer
// Outputs:
//   p_buffer: value set to newly allocated buffer that has been initialized with file contents
//   returned value: length of allocated buffer
// Notes:
//   Must be run on little-endian machine.
//   Implements NumPy .npy format version 1.0:
//       numpy.org/doc/stable/reference/generated/numpy.lib.format.html#module-numpy.lib.format
//   A .npy file is simply a binary dump of a NumPy array prepended with an ASCII string that can be
//       interpreted as a Python dictionary literal with certain entries.
std::size_t OutputWriter::GenerateNpyFromDoubleArray(const Array<double> &array, uint8_t **p_buffer)
{
  // Prepare buffer
  const std::size_t header_length = 128;
  std::size_t data_length = array.GetNumBytes();
  std::size_t buffer_length = header_length + data_length;
  *p_buffer = new uint8_t[buffer_length];

  // Write magic string and version number to buffer
  std::size_t length = 0;
  const char *magic_string = "\x93NUMPY";
  std::memcpy(*p_buffer + length, magic_string, 6);
  length += 6;
  const char *version = "\x01\x00";
  std::memcpy(*p_buffer + length, version, 2);
  length += 2;

  // Write length of header proper to buffer
  uint16_t header_len = static_cast<uint16_t>(header_length - length - 2);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&header_len), 2);
  length += 2;

  // Write header proper to buffer
  char *buffer_address = reinterpret_cast<char *>(*p_buffer + length);
  int num_written = -1;
  if (array.n2 == 1 and array.n3 == 1 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d,)}", array.n1);
  else if (array.n3 == 1 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d)}", array.n2, array.n1);
  else if (array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d, %d)}", array.n3, array.n2,
        array.n1);
  else if (array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d, %d, %d)}", array.n4, array.n3,
        array.n2, array.n1);
  else
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d, %d, %d, %d)}", array.n5,
        array.n4, array.n3, array.n2, array.n1);
  if (num_written < 0 or num_written > static_cast<int>(header_length - length - 1))
    throw BlacklightException("Error converting data to .npy format.");
  std::size_t num_spaces =
      static_cast<std::size_t>(static_cast<int>(header_length - length - 1) - num_written);
  if (num_spaces > 0)
    std::memset(buffer_address + num_written, ' ', num_spaces);
  (*p_buffer)[header_length-1] = '\n';
  length = header_length;

  // Write array to buffer
  const uint8_t *data_pointer = reinterpret_cast<const uint8_t *>(array.data);
  std::memcpy(*p_buffer + length, data_pointer, data_length);
  length += data_length;
  return buffer_length;
}
