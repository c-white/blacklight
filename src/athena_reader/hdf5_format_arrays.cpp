// Blacklight Athena++ reader - HDF5 interface for reading arrays

// C++ headers
#include <cstring>  // memcpy
#include <string>   // string

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "athena_reader.hpp"
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function to read integer array dataset from HDF5 file by name
// Inputs:
//   name: name of dataset
// Outputs:
//   int_array: array allocated and set (indirectly)
// Notes:
//   Changes stream pointer.
void AthenaReader::ReadHDF5IntArray(const char *name, Array<int> &int_array)
{
  // Locate header
  unsigned long int header_address = ReadHDF5DatasetHeaderAddress(name);

  // Read header
  unsigned char *datatype_raw, *dataspace_raw, *data_raw;
  ReadHDF5DataObjectHeader(header_address, &datatype_raw, &dataspace_raw, &data_raw);

  // Set array
  SetHDF5IntArray(datatype_raw, dataspace_raw, data_raw, int_array);
  delete[] datatype_raw;
  delete[] dataspace_raw;
  delete[] data_raw;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read float array dataset from HDF5 file by name
// Inputs:
//   name: name of dataset
// Outputs:
//   float_array: array allocated and set (indirectly)
// Notes:
//   Changes stream pointer.
void AthenaReader::ReadHDF5FloatArray(const char *name, Array<float> &float_array)
{
  // Locate header
  unsigned long int header_address = ReadHDF5DatasetHeaderAddress(name);

  // Read header
  unsigned char *datatype_raw, *dataspace_raw, *data_raw;
  ReadHDF5DataObjectHeader(header_address, &datatype_raw, &dataspace_raw, &data_raw);

  // Set array
  SetHDF5FloatArray(datatype_raw, dataspace_raw, data_raw, float_array);
  delete[] datatype_raw;
  delete[] dataspace_raw;
  delete[] data_raw;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to initialize string array dataset from HDF5
// Inputs:
//   datatype_raw: raw datatype description
//   dataspace_raw: raw dataspace description
//   data_raw: raw data
//   allocate: flag indicating new memory should be allocated
// Outputs:
//   *p_string_array: array allocated and set
//   *p_array_length: number of allocated members
// Notes:
//   Must have datatype version 1.
//   Must be a 1D array.
void AthenaReader::SetHDF5StringArray(const unsigned char *datatype_raw,
    const unsigned char *dataspace_raw, const unsigned char *data_raw, bool allocate,
    std::string **p_string_array, int *p_array_length)
{
  // Check datatype version and class
  int offset = 0;
  unsigned char version_class = datatype_raw[offset++];
  if (version_class >> 4 != 1)
    throw BlacklightException("Unexpected HDF5 datatype version.");
  if ((version_class & 0b00001111) != 3)
    throw BlacklightException("Unexpected HDF5 datatype class.");

  // Read datatype metadata
  unsigned char class_1 = datatype_raw[offset++];
  offset += 2;

  // Read data size
  unsigned int size;
  std::memcpy(&size, datatype_raw + offset, 4);
  offset += 4;

  // Check character set
  if (class_1 >> 4 != 0)
    throw BlacklightException("Unexpected HDF5 string encoding.");

  // Read and check dimensions
  unsigned long int *dims;
  int num_dims;
  ReadHDF5DataspaceDims(dataspace_raw, &dims, &num_dims);
  if (num_dims != 1)
    throw BlacklightException("Unexpected HDF5 string array size.");
  if (allocate)
    *p_array_length = static_cast<int>(dims[0]);
  else if (*p_array_length != static_cast<int>(dims[0]))
    throw BlacklightException("Array mismatch upon subsequent read.");

  // Allocate and initialize array
  if (allocate)
    *p_string_array = new std::string[*p_array_length];
  char *buffer = new char[size];
  for (int n = 0; n < *p_array_length; n++)
  {
    std::memcpy(buffer, data_raw + size * static_cast<unsigned int>(n), size);
    (*p_string_array)[n] = buffer;
  }
  delete[] buffer;

  // Free dimensions
  delete[] dims;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to initialize integer array dataset from HDF5
// Inputs:
//   datatype_raw: raw datatype description
//   dataspace_raw: raw dataspace description
//   data_raw: raw data
// Outputs:
//   int_array: array allocated (if not already allocated) and set
// Notes:
//   Must have datatype version 1.
//   Must be 4 or 8 bytes.
//   8-byte integers will be truncated to 4 bytes.
//   Must have trivial padding.
//   Must have no offset.
//   Must be run on little-endian machine.
void AthenaReader::SetHDF5IntArray(const unsigned char *datatype_raw,
    const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<int> &int_array)
{
  // Check datatype version and class
  int offset = 0;
  unsigned char version_class = datatype_raw[offset++];
  if (version_class >> 4 != 1)
    throw BlacklightException("Unexpected HDF5 datatype version.");
  if ((version_class & 0b00001111) != 0)
    throw BlacklightException("Unexpected HDF5 datatype class.");

  // Read datatype metadata
  unsigned char class_1 = datatype_raw[offset++];
  offset += 2;

  // Read data size
  unsigned int size;
  std::memcpy(&size, datatype_raw + offset, 4);
  offset += 4;
  if (size != 4 and size != 8)
    throw BlacklightException("Unexpected int size.");

  // Read and check properties
  bool rev_endian = class_1 & 0b00000001;
  if (class_1 & 0b00000010 or class_1 & 0b00000100)
    throw BlacklightException("Unexpected HDF5 fixed-point padding.");
  bool signed_val = class_1 & 0b00001000;
  unsigned short int bit_offset, bit_precision;
  std::memcpy(&bit_offset, datatype_raw + offset, 2);
  offset += 2;
  std::memcpy(&bit_precision, datatype_raw + offset, 2);
  offset += 2;
  if (bit_offset != 0 or bit_precision != 8 * size)
    throw BlacklightException("Unexpected HDF5 fixed-point bit layout.");

  // Read dimensions
  unsigned long int *dims;
  int num_dims;
  ReadHDF5DataspaceDims(dataspace_raw, &dims, &num_dims);
  unsigned int num_elements = 1;
  for (int n = 0; n < num_dims; n++)
    num_elements *= static_cast<unsigned int>(dims[n]);

  // Allocate array
  if (num_dims == 1)
  {
    if (not int_array.allocated)
      int_array.Allocate(static_cast<int>(dims[0]));
    else if (static_cast<unsigned long int>(int_array.n1) != dims[0]
        or static_cast<unsigned int>(int_array.n_tot) != num_elements)
      throw BlacklightException("Array mismatch upon subsequent read.");
  }
  else if (num_dims == 2)
  {
    if (not int_array.allocated)
      int_array.Allocate(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
    else if (static_cast<unsigned long int>(int_array.n2) != dims[0]
        or static_cast<unsigned long int>(int_array.n1) != dims[1]
        or static_cast<unsigned int>(int_array.n_tot) != num_elements)
      throw BlacklightException("Array mismatch upon subsequent read.");
  }
  else
    throw BlacklightException("Unexpected HDF5 fixed-point array size.");
  delete[] dims;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate buffer
    char *buffer = new char[size];

    // Initialize array
    #pragma omp for schedule(static)
    for (unsigned int n = 0; n < num_elements; n++)
    {
      if (rev_endian)
        for (unsigned int m = 0; m < size; m++)
          std::memcpy(buffer + size - 1 - m, data_raw + n * size + m, 1);
      else
        std::memcpy(buffer, data_raw + n * size, 4);
      if (signed_val)
        int_array(n) = *reinterpret_cast<int *>(buffer);
      else
        int_array(n) = static_cast<int>(*reinterpret_cast<unsigned int *>(buffer));
    }

    // Free buffer
    delete[] buffer;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to initialize 4-byte floating point array dataset from HDF5
// Inputs:
//   datatype_raw: raw datatype description
//   dataspace_raw: raw dataspace description
//   data_raw: raw data
// Outputs:
//   float_array: array allocated (if not already allocated) and set
// Notes:
//   Must have datatype version 1.
//   Must be standard 4-byte floats.
//   Must be run on little-endian machine.
void AthenaReader::SetHDF5FloatArray(const unsigned char *datatype_raw,
    const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array)
{
  // Check datatype version and class
  int offset = 0;
  unsigned char version_class = datatype_raw[offset++];
  if (version_class >> 4 != 1)
    throw BlacklightException("Unexpected HDF5 datatype version.");
  if ((version_class & 0b00001111) != 1)
    throw BlacklightException("Unexpected HDF5 datatype class.");

  // Read datatype metadata
  unsigned char class_1 = datatype_raw[offset++];
  unsigned char class_2 = datatype_raw[offset++];
  offset++;

  // Read data size
  unsigned int size;
  std::memcpy(&size, datatype_raw + offset, 4);
  offset += 4;
  if (size != 4)
    throw BlacklightException("Unexpected float size.");

  // Check properties
  bool rev_endian = class_1 & 0b00000001;
  if (class_1 & 0b01000000)
    throw BlacklightException("Unexpected HDF5 floating-point byte order.");
  if (class_1 & 0b00001110)
    throw BlacklightException("Unexpected HDF5 floating-point padding.");
  if ((class_1 & 0b00110000) != 0b00100000)
    throw BlacklightException("Unexpected HDF5 floating-point mantissa normalization.");
  if (class_2 != 31)
    throw BlacklightException("Unexpected HDF5 floating-point sign location.");
  unsigned short int bit_offset, bit_precision;
  unsigned char exp_loc, exp_size, man_loc, man_size;
  unsigned int exp_bias;
  std::memcpy(&bit_offset, datatype_raw + offset, 2);
  offset += 2;
  std::memcpy(&bit_precision, datatype_raw + offset, 2);
  offset += 2;
  std::memcpy(&exp_loc, datatype_raw + offset++, 1);
  std::memcpy(&exp_size, datatype_raw + offset++, 1);
  std::memcpy(&man_loc, datatype_raw + offset++, 1);
  std::memcpy(&man_size, datatype_raw + offset++, 1);
  std::memcpy(&exp_bias, datatype_raw + offset, 4);
  offset += 4;
  if (bit_offset != 0 or bit_precision != 32 or exp_loc != 23 or exp_size != 8 or man_loc != 0
      or man_size != 23 or exp_bias != 127)
    throw BlacklightException("Unexpected HDF5 floating-point bit layout.");

  // Read dimensions
  unsigned long int *dims;
  int num_dims;
  ReadHDF5DataspaceDims(dataspace_raw, &dims, &num_dims);
  unsigned int num_elements = 1;
  for (int n = 0; n < num_dims; n++)
    num_elements *= static_cast<unsigned int>(dims[n]);

  // Allocate array
  if (num_dims == 0)
  {
    if (not float_array.allocated)
      float_array.Allocate(1);
    else if (static_cast<unsigned int>(float_array.n_tot) != num_elements)
      throw BlacklightException("Array mismatch upon subsequent read.");
  }
  else if (num_dims == 1)
  {
    if (not float_array.allocated)
      float_array.Allocate(static_cast<int>(dims[0]));
    else if (static_cast<unsigned long int>(float_array.n1) != dims[0]
        or static_cast<unsigned int>(float_array.n_tot) != num_elements)
      throw BlacklightException("Array mismatch upon subsequent read.");
  }
  else if (num_dims == 2)
  {
    if (not float_array.allocated)
      float_array.Allocate(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
    else if (static_cast<unsigned long int>(float_array.n2) != dims[0]
        or static_cast<unsigned long int>(float_array.n1) != dims[1]
        or static_cast<unsigned int>(float_array.n_tot) != num_elements)
      throw BlacklightException("Array mismatch upon subsequent read.");
  }
  else if (num_dims == 5)
  {
    if (not float_array.allocated)
      float_array.Allocate(static_cast<int>(dims[0]), static_cast<int>(dims[1]),
          static_cast<int>(dims[2]), static_cast<int>(dims[3]), static_cast<int>(dims[4]));
    else if (static_cast<unsigned long int>(float_array.n5) != dims[0]
        or static_cast<unsigned long int>(float_array.n4) != dims[1]
        or static_cast<unsigned long int>(float_array.n3) != dims[2]
        or static_cast<unsigned long int>(float_array.n2) != dims[3]
        or static_cast<unsigned long int>(float_array.n1) != dims[4]
        or static_cast<unsigned int>(float_array.n_tot) != num_elements)
      throw BlacklightException("Array mismatch upon subsequent read.");
  }
  else
    throw BlacklightException("Unexpected HDF5 floating-point array size.");
  delete[] dims;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate buffer
    char *buffer = new char[size];

    // Initialize array
    #pragma omp for schedule(static)
    for (unsigned int n = 0; n < num_elements; n++)
    {
      if (rev_endian)
        for (unsigned int m = 0; m < size; m++)
          std::memcpy(buffer + size - 1 - m, data_raw + n * size + m, 1);
      else
        std::memcpy(buffer, data_raw + n * size, size);
      float_array(n) = *reinterpret_cast<float *>(buffer);
    }

    // Free buffer
    delete[] buffer;
  }
  return;
}
