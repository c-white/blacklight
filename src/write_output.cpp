// Blacklight output writer

// C++ headers
#include <cstdint>   // uint8_t, uint16_t, uint32_t
#include <cstdio>    // size_t, snprintf
#include <cstring>   // memcpy, memset
#include <ctime>     // localtime, time, time_t, tm
#include <fstream>   // ofstream
#include <ios>       // ios_base, streamsize
#include <optional>  // optional

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

//--------------------------------------------------------------------------------------------------

// Output writer for a raw binary dump without metadata
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
//   Output file endianness will match that of machine writing file.
void OutputWriter::WriteRaw()
{
  const char *data_pointer = reinterpret_cast<const char *>(image.data);
  std::streamsize data_size = static_cast<std::streamsize>(image.GetNumBytes());
  p_output_stream->write(data_pointer, data_size);
  return;
}

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

//--------------------------------------------------------------------------------------------------

// Function for populating a buffer with a ZIP local file header
// Inputs:
//   record: buffer containing entire record to be written to ZIP archive
//   record_length: size of record
//   record_name: name of record to encode in header
// Outputs:
//   p_buffer: value set to newly allocated buffer that has been initialized with header contents
//   returned value: length of allocated buffer
// Notes:
//   Must be run on little-endian machine.
//   Name in file will be given record_name + ".npy".
std::size_t OutputWriter::GenerateZIPLocalFileHeader(const uint8_t *record,
    std::size_t record_length, const char *record_name, uint8_t **p_buffer)
{
  // Prepare buffer
  std::size_t name_length = std::strlen(record_name) + 4;
  std::size_t buffer_length = 30 + name_length;
  *p_buffer = new uint8_t[buffer_length];

  // Write signature to buffer
  std::size_t length = 0;
  const char *signature = "\x50\x4b\x03\x04";
  std::memcpy(*p_buffer, signature, 4);
  length += 4;

  // Write version needed to buffer
  const uint8_t version_needed_compatibility = 0;
  const uint8_t version_needed_major = 2;
  const uint8_t version_needed_minor = 0;
  const uint8_t version_needed_both = 10 * version_needed_major + version_needed_minor;
  std::memcpy(*p_buffer + length, &version_needed_both, 1);
  length += 1;
  std::memcpy(*p_buffer + length, &version_needed_compatibility, 1);
  length += 1;

  // Write bit flags to buffer
  const uint16_t bit_flags = 0;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&bit_flags), 2);
  length += 2;

  // Write compression method to buffer
  const uint16_t compression = 0;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&compression), 2);
  length += 2;

  // Write current date and time to buffer
  std::time_t current_time = std::time(nullptr);
  std::tm *current_local_time = std::localtime(&current_time);
  uint16_t time = static_cast<uint16_t>(current_local_time->tm_hour << 11);
  time |= static_cast<uint16_t>((current_local_time->tm_min & 0b0000000000011111) << 5);
  time |= static_cast<uint16_t>(current_local_time->tm_sec / 2 & 0b0000000000011111);
  uint16_t date = static_cast<uint16_t>((current_local_time->tm_year - 80) << 9);
  date |= static_cast<uint16_t>(((current_local_time->tm_mon + 1) & 0b0000000000001111) << 5);
  date |= static_cast<uint16_t>(current_local_time->tm_mday & 0b0000000000011111);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&time), 2);
  length += 2;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&date), 2);
  length += 2;

  // Write CRC-32 to buffer
  uint32_t crc_32 = CalculateCRC32(record, record_length);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&crc_32), 4);
  length += 4;

  // Write compressed and uncompressed sizes to buffer
  if (record_length > UINT32_MAX)
    throw BlacklightException("Array and metadata too large for ZIP record.");
  uint32_t record_length_32 = static_cast<uint32_t>(record_length);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&record_length_32), 4);
  length += 4;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&record_length_32), 4);
  length += 4;

  // Write name length and extra field length to buffer
  if (name_length > UINT16_MAX)
    throw BlacklightException("Array name too long for ZIP format.");
  uint16_t name_length_16 = static_cast<uint16_t>(name_length);
  const uint16_t extra_field_length = 0;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&name_length_16), 2);
  length += 2;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&extra_field_length), 2);
  length += 2;

  // Write name to buffer
  const char *record_name_extension = ".npy";
  std::memcpy(*p_buffer + length, record_name, name_length - 4);
  length += name_length - 4;
  std::memcpy(*p_buffer + length, record_name_extension, 4);
  length += 4;
  return buffer_length;
}

//--------------------------------------------------------------------------------------------------

// Function for populating a buffer with a ZIP central directory header
// Inputs:
//   local_header: buffer containing local file header
//   local_header_offset: offset of start of local file header from start of file
// Outputs:
//   p_buffer: value set to newly allocated buffer that has been initialized with header contents
//   returned value: length of allocated buffer
// Notes:
//   Must be run on little-endian machine.
//   File attributes are not fully specified by ZIP standard. The values written match what NumPy
//       writes on at least one system.
std::size_t OutputWriter::GenerateZIPCentralDirectoryHeader(const uint8_t *local_header,
    std::size_t local_header_offset, uint8_t **p_buffer)
{
  // Extract length of file name
  uint16_t name_length = *reinterpret_cast<const uint16_t *>(local_header + 26);

  // Prepare buffer
  std::size_t buffer_length = static_cast<std::size_t>(46 + name_length);
  *p_buffer = new uint8_t[buffer_length];

  // Write signature to buffer
  std::size_t length = 0;
  const char *signature = "\x50\x4b\x01\x02";
  std::memcpy(*p_buffer, signature, 4);
  length += 4;

  // Write version made by to buffer
  const uint8_t version_made_compatibility = 3;
  const uint8_t version_made_major = 2;
  const uint8_t version_made_minor = 0;
  const uint8_t version_made_both = 10 * version_made_major + version_made_minor;
  std::memcpy(*p_buffer + length, &version_made_both, 1);
  length += 1;
  std::memcpy(*p_buffer + length, &version_made_compatibility, 1);
  length += 1;

  // Copy appropriate fields from local file header to buffer
  std::memcpy(*p_buffer + length, local_header + 4, 26);
  length += 26;

  // Write file comment length to buffer
  const uint16_t comment_length = 0;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&comment_length), 2);
  length += 2;

  // Write disk number start to buffer
  const uint16_t disk_number = 0;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&disk_number), 2);
  length += 2;

  // Write file attributes to buffer
  const uint16_t internal_attributes = 0;
  const uint32_t external_attributes = 0b10000001100000000000000000000000;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&internal_attributes), 2);
  length += 2;
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&external_attributes), 4);
  length += 4;

  // Write local header offset to buffer
  uint32_t local_header_offset_32 = static_cast<uint32_t>(local_header_offset);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&local_header_offset_32), 4);
  length += 4;

  // Write name to buffer
  std::memcpy(*p_buffer + length, local_header + 30, name_length);
  length += name_length;
  return buffer_length;
}

//--------------------------------------------------------------------------------------------------

// Function for populating a buffer with a ZIP end of central directory record
// Inputs:
//   central_directory_offset: offset of first central directory file header from start of file
//   central_directory_length: size of central directory, excluding this record
//   num_central_directory_entries: number of entries in central directory
// Outputs:
//   p_buffer: value set to newly allocated buffer that has been initialized with header contents
//   returned value: length of allocated buffer
// Notes:
//   Must be run on little-endian machine.
std::size_t OutputWriter::GenerateZIPEndOfCentralDirectoryRecord(
    std::size_t central_directory_offset, std::size_t central_directory_length,
    int num_central_directory_entries, uint8_t **p_buffer)
{
  // Prepare buffer
  const std::size_t buffer_length = 22;
  *p_buffer = new uint8_t[buffer_length];

  // Write signature to buffer
  std::size_t length = 0;
  const char *signature = "\x50\x4b\x05\x06";
  std::memcpy(*p_buffer, signature, 4);
  length += 4;

  // Write disk numbers to buffer
  const uint16_t disk_number = 0;
  std::memcpy(*p_buffer + length, &disk_number, 2);
  length += 2;
  std::memcpy(*p_buffer + length, &disk_number, 2);
  length += 2;

  // Write number of entries to buffer
  uint16_t num_central_directory_entries_16 = static_cast<uint16_t>(num_central_directory_entries);
  std::memcpy(*p_buffer + length, &num_central_directory_entries_16, 2);
  length += 2;
  std::memcpy(*p_buffer + length, &num_central_directory_entries_16, 2);
  length += 2;

  // Write central directory size to buffer
  if (central_directory_length > UINT32_MAX)
    throw BlacklightException("Central directory too large for ZIP format.");
  uint32_t central_directory_length_32 = static_cast<uint32_t>(central_directory_length);
  std::memcpy(*p_buffer + length, &central_directory_length_32, 4);
  length += 4;

  // Write central directory offset to buffer
  if (central_directory_offset > UINT32_MAX)
    throw BlacklightException("File contents too large for ZIP format.");
  uint32_t central_directory_offset_32 = static_cast<uint32_t>(central_directory_offset);
  std::memcpy(*p_buffer + length, &central_directory_offset_32, 4);
  length += 4;

  // Write file comment length to buffer
  const uint16_t file_comment_length = 0;
  std::memcpy(*p_buffer + length, &file_comment_length, 2);
  length += 2;
  return buffer_length;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating cyclic redundancy check on a message
// Inputs:
//   message: message for which value should be calculated
//   message_length: number of bytes in message
// Outputs:
//   returned value: CRC-32 of message
// Notes:
//   Calculates CRC-32 as defined for ZIP format, effectively doing the following:
//     Takes as input a message M understood to be first byte first, most significant bit first
//         within each byte.
//     Uses the degree-32 polynomial P(x) = x^32 + x^26 + x^23 + x^22 + x^16 + x^12 + x^11 + x^10
//         + x^8 + x^7 + x^5 + x^4 + x^2 + x^1 + x^0. Omitting the leading term, this can be encoded
//         in binary as 0b00000100110000010001110110110111 (leading coefficient (possibly 0) first,
//         most significant bit first, no notion of bytes). This in turn can be encoded in
//         hexadecimal as 0x04C11DB7 (most significant value first, no notion of bytes). This
//         polynomial is not explicitly given in the ZIP standard.
//     Calculates M' as M with each byte bit-reversed (first byte still first, least significant bit
//         first within each byte).
//     Appends 32 0 bits to the end of M'.
//     Inverts the first 32 bits of M' for preconditioning.
//     Interprets M' as a polynomial of degree one less than its length, with the first bit
//         corresponding to the leading coefficient (possibly 0) and with no notion of bytes.
//     Calculates the degree-31 polynomial remainder R = M' % P over the field with 2 elements.
//     Interprets R as 32 bits, with the leading coefficient (possibly 0) first and with no notion
//         of bytes.
//     Calculates R' as R bit-reversed. Unlike with calculating M' from M, here there is no notion
//         of bytes.
//     Inverts R' for postconditioning.
//     Calculates R'' as R' byte-reversed. That is, if R is taken to be least significant bit first
//         (like M'), then R' can be considered most significant bit first, or big-endian with most
//         significant bit first within each byte. Then R'' has the same value as R' except in a
//         little-endian form with most significant bit first within each byte.
//     Returns R'', the little-endian, 4-byte int version of the remainder.
//   As an example, if M consists of 32 0 bits, R' has the value 0x2144DF1C = 558161692, while R''
//       has the bit pattern 0x1CDF4421 = 0b00011100110111110100010000100001.
//   The ZIP standard refers to a "magic number" 0xDEBB20E3. This can be obtained as follows for any
//       M:
//     Calculate M' by bit-reversing M as above.
//     Append inverted R (that is, invert all the bits of R as per the postconditioning above but do
//         not reverse them) to M'.
//     Append 32 0 bits to the end of M' as above.
//     Invert the first 32 bits of M' as above.
//     Calculate a new remainder S = M' % P as above.
//     Calculate S' as S bit-reversed with no notion of bytes.
//     S' will have the bit pattern 0xDEBB20E3 = 0b11011110101110110010000011100011.
//   Must be run on little-endian machine.
//   Does not handle messages less than 4 bytes long.
uint32_t OutputWriter::CalculateCRC32(const uint8_t *message, std::size_t message_length)
{
  // Check message length
  if (message_length < 4)
    throw BlacklightException("CRC-32 could not be calculated.");

  // Calculate lookup table for reversing bits within a byte
  uint8_t bit_reverse_table[256];
  for (int input = 0; input < 256; input++)
  {
    uint8_t byte = static_cast<uint8_t>(input);
    byte = static_cast<uint8_t>((byte & 0b11110000) >> 4 | (byte & 0b00001111) << 4);
    byte = static_cast<uint8_t>((byte & 0b11001100) >> 2 | (byte & 0b00110011) << 2);
    byte = static_cast<uint8_t>((byte & 0b10101010) >> 1 | (byte & 0b01010101) << 1);
    bit_reverse_table[input] = byte;
  }

  // Calculate lookup table for dividing by polynomial
  const uint64_t polynomial = 0x04C11DB700000000ULL;
  uint32_t division_table[256] = {};
  for (int input = 0; input < 256; input++)
  {
    uint64_t scratch = static_cast<uint64_t>(input) << 56;
    if (scratch & 0x8000000000000000ULL)
      scratch ^= polynomial >> 1;
    if (scratch & 0x4000000000000000ULL)
      scratch ^= polynomial >> 2;
    if (scratch & 0x2000000000000000ULL)
      scratch ^= polynomial >> 3;
    if (scratch & 0x1000000000000000ULL)
      scratch ^= polynomial >> 4;
    if (scratch & 0x0800000000000000ULL)
      scratch ^= polynomial >> 5;
    if (scratch & 0x0400000000000000ULL)
      scratch ^= polynomial >> 6;
    if (scratch & 0x0200000000000000ULL)
      scratch ^= polynomial >> 7;
    if (scratch & 0x0100000000000000ULL)
      scratch ^= polynomial >> 8;
    uint32_t output = static_cast<uint32_t>(scratch >> 24 & 0x00000000FFFFFFFFULL);
    output = (output & 0xFFFF0000) >> 16 | (output & 0x0000FFFF) << 16;
    output = (output & 0xFF00FF00) >> 8 | (output & 0x00FF00FF) << 8;
    division_table[input] = output;
  }

  // Perform division on preconditioned, bit-reversed message appended with 0's
  uint8_t leading_byte;
  uint32_t remainder = 0xFFFFFFFF;
  uint8_t *remainder_bytes = reinterpret_cast<uint8_t *>(&remainder);
  for (int b = 0; b < 4; b++)
    remainder_bytes[b] ^= bit_reverse_table[message[b]];
  for (std::size_t n = 4; n < message_length; n++)
  {
    leading_byte = remainder_bytes[0];
    for (int b = 0; b < 3; b++)
      remainder_bytes[b] = remainder_bytes[b+1];
    remainder_bytes[3] = bit_reverse_table[message[n]];
    remainder ^= division_table[leading_byte];
  }
  for (int n = 0; n < 4; n++)
  {
    leading_byte = remainder_bytes[0];
    for (int b = 0; b < 3; b++)
      remainder_bytes[b] = remainder_bytes[b+1];
    remainder_bytes[3] = 0;
    remainder ^= division_table[leading_byte];
  }

  // Modify remainder to account for bit reversing and postconditioning
  for (int n = 0; n < 4; n++)
    remainder_bytes[n] = bit_reverse_table[remainder_bytes[n]];
  remainder ^= 0xFFFFFFFF;
  return remainder;
}
