// Blacklight output writer - ZIP output format

// C++ headers
#include <cstdint>  // uint8_t, uint16_t, uint32_t
#include <cstdio>   // size_t
#include <cstring>  // memcpy
#include <ctime>    // localtime, time, time_t, tm

// Blacklight headers
#include "output_writer.hpp"
#include "../utils/exceptions.hpp"  // BlacklightException

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
