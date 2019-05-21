// Ray Trace Athena++ reader

// C++ headers
#include <cstring>  // memcpy, size_t
#include <fstream>  // ifstream
#include <ios>      // streamoff
#include <string>   // string

// Ray Trace headers
#include "read_athena.hpp"
#include "read_input.hpp"   // input_reader
#include "exceptions.hpp"   // ray_trace_exception

//--------------------------------------------------------------------------------------------------

// Athena++ reader constructor
// Inputs:
//   data_file: name of input file
// Notes:
//   Opens stream for reading.
athena_reader::athena_reader(const std::string data_file)
  : data_stream(data_file)
{
  // Check that file is open
  if (not data_stream.is_open())
    throw ray_trace_exception("Error: Could not open data file.\n");
}

//--------------------------------------------------------------------------------------------------

// Athena++ reader destructor
athena_reader::~athena_reader() {}

//--------------------------------------------------------------------------------------------------

// Athena++ reader read and initialize function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Initializes all member objects.
void athena_reader::read()
{
  // Read superblock
  read_hdf5_superblock();

  // Read necessary file attributes
  read_hdf5_root_object_header();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 superblock
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Indirectly sets root_object_header_address, btree_address, and root_name_heap_address.
//   Changes stream pointer.
//   Does not allow for superblock to be anywhere but beginning of file.
//   Must have superblock version 0.
//   Must have size of offsets 8.
//   Must have size of lengths 8.
void athena_reader::read_hdf5_superblock()
{
  // Check format signature
  data_stream.seekg(0);
  const unsigned char expected_signature[] = {0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a};
  for (int n = 0; n < 8; n++)
    if (data_stream.get() != expected_signature[n])
      throw ray_trace_exception("Error: Unexpected HDF5 format signature.\n");

  // Check superblock version
  if (data_stream.get() != 0)
    throw ray_trace_exception("Error: Unexpected HDF5 superblock version.\n");

  // Check other version numbers
  if (data_stream.get() != 0)
    throw ray_trace_exception("Error: Unexpected HDF5 file free space storage version.\n");
  if (data_stream.get() != 0)
    throw ray_trace_exception("Error: Unexpected HDF5 root group symbol table entry version.\n");
  data_stream.ignore(1);
  if (data_stream.get() != 0)
    throw ray_trace_exception("Error: Unexpected HDF5 shared header message format version.\n");

  // Check sizes
  if (data_stream.get() != 8)
    throw ray_trace_exception("Error: Unexpected HDF5 size of offsets.\n");
  if (data_stream.get() != 8)
    throw ray_trace_exception("Error: Unexpected HDF5 size of lengths.\n");
  data_stream.ignore(1);

  // Skip checking tree parameters and consistency flags
  data_stream.ignore(2 * 2 + 4);

  // Skip checking addresses
  data_stream.ignore(4 * 8);

  // Read root group symbol table entry
  read_root_group_symbol_table_entry();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root group symbol table entry
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets root_object_header_address, btree_address, and root_name_heap_address.
//   Assumes stream pointer is already set.
//   Must have size of offsets 8.
//   Must be run on little-endian machine
void athena_reader::read_root_group_symbol_table_entry()
{
  // Skip reading link name offset
  data_stream.ignore(8);

  // Read object header address
  data_stream.read(reinterpret_cast<char *>(&root_object_header_address), 8);

  // Check cache type
  unsigned int cache_type;
  data_stream.read(reinterpret_cast<char *>(&cache_type), 4);
  if (cache_type != 1)
    throw ray_trace_exception("Error: Unexpected HDF5 root group symbol table entry cache type.\n");

  // Read B-tree and name heap addresses
  data_stream.read(reinterpret_cast<char *>(&btree_address), 8);
  data_stream.read(reinterpret_cast<char *>(&root_name_heap_address), 8);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root object header
// Inputs:
//   data_stream: input file stream
//   root_object_header_address: offset into file of root header
// Outputs:
//   *p_dataset_names: newly allocated array of names of cell value datasets
//   *p_variable_names: newly allocated array of names of variables in cell value datasets
//   *p_num_variables: newly allocated array of numbers of variables in each dataset
// Notes:
//   Allocates and sets dataset_names, variable_names, nums_variables.
//   Sets num_dataset_names and num_variable_names.
//   Assumes root_object_header_address set.
//   Changes stream pointer.
//   Must have object header version 1.
//   Must not have shared header messages.
//   Must have attribute message version 1.
void athena_reader::read_hdf5_root_object_header()
{
  // Check object header version
  data_stream.seekg(static_cast<std::streamoff>(root_object_header_address));
  if (data_stream.get() != 1)
    throw ray_trace_exception("Error: Unexpected HDF5 object header version.\n");
  data_stream.ignore(1);

  // Read number of header messages
  unsigned short int num_messages;
  data_stream.read(reinterpret_cast<char *>(&num_messages), 2);

  // Skip reading object reference count and object header size
  data_stream.ignore(8);

  // Align to 8 bytes within header (location of padding not documented)
  data_stream.ignore(4);

  // Go through messages
  bool dataset_names_found = false;
  bool variable_names_found = false;
  bool num_variables_found = false;
  for (int n = 0; n < num_messages; n++)
  {
    // Read message type and size
    unsigned short int message_type, message_size;
    data_stream.read(reinterpret_cast<char *>(&message_type), 2);
    data_stream.read(reinterpret_cast<char *>(&message_size), 2);

    // Check message flags
    unsigned char message_flags;
    data_stream.read(reinterpret_cast<char *>(&message_flags), 1);
    data_stream.ignore(3);
    if (message_flags & 0b00000010)
      throw ray_trace_exception("Error: Unexpected HDF5 header message flag.\n");

    // Read message data
    unsigned char *message_data = new unsigned char[message_size];
    data_stream.read(reinterpret_cast<char *>(message_data), message_size);

    // Follow any continuation messages
    if (message_type == 16)
    {
      unsigned long int new_offset;
      std::memcpy(&new_offset, message_data, 8);
      data_stream.seekg(static_cast<std::streamoff>(new_offset));
      continue;

    // Inspect any attribute messages
    } else if (message_type == 12) {

      // Check attribute message version
      int offset = 0;
      if (message_data[0] != 1)
        throw ray_trace_exception("Error: Unexpected HDF5 attribute message version.\n");
      offset += 2;

      // Read attribute message metadata
      unsigned short int name_size, datatype_size, dataspace_size;
      std::memcpy(&name_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&datatype_size, message_data + offset, 2);
      offset += 2;
      std::memcpy(&dataspace_size, message_data + offset, 2);
      offset += 2;
      unsigned short int name_size_pad = static_cast<unsigned short int>((8 - name_size % 8) % 8);
      unsigned short int datatype_size_pad =
          static_cast<unsigned short int>((8 - datatype_size % 8) % 8);
      unsigned short int dataspace_size_pad =
          static_cast<unsigned short int>((8 - dataspace_size % 8) % 8);

      // Read attribute message data
      unsigned char *name_raw = new unsigned char[name_size];
      unsigned char *datatype_raw = new unsigned char[datatype_size];
      unsigned char *dataspace_raw = new unsigned char[dataspace_size];
      std::memcpy(name_raw, message_data + offset, name_size);
      offset += name_size + name_size_pad;
      std::memcpy(datatype_raw, message_data + offset, datatype_size);
      offset += datatype_size + datatype_size_pad;
      std::memcpy(dataspace_raw, message_data + offset, dataspace_size);
      offset += dataspace_size + dataspace_size_pad;
      std::string name(reinterpret_cast<char *>(name_raw),
          static_cast<std::string::size_type>(static_cast<int>(name_size) - 1));

      // Read and set desired attributes
      if (name == "DatasetNames")
      {
        dataset_names_found = true;
        set_hdf5_string_array(datatype_raw, dataspace_raw, message_data + offset, &dataset_names,
            &num_dataset_names);
      } else if (name == "VariableNames") {
        variable_names_found = true;
      } else if (name == "NumVariables") {
        num_variables_found = true;
      }

      // Free raw buffers
      delete[] name_raw;
      delete[] datatype_raw;
      delete[] dataspace_raw;
    }

    // Free raw buffer
    delete[] message_data;

    // Break when required information found
    if (dataset_names_found and variable_names_found and num_variables_found)
      break;
  }

  // Check that appropriate messages were found
  if (not (dataset_names_found and variable_names_found and num_variables_found))
    throw ray_trace_exception("Error: Could not find needed file-level attributes.\n");
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to initialize string array dataset from HDF5
// Inputs:
//   datatype_raw: raw datatype description
//   dataspace_raw: raw dataspace description
//   data_raw: raw data
// Outputs:
//   *p_string_array: array allocated and initialized
//   *p_array_length: number of allocated members
// Notes:
//   Must have datatype version 1.
//   Must be a 1D array.

void athena_reader::set_hdf5_string_array(const unsigned char *datatype_raw,
    const unsigned char *dataspace_raw, const unsigned char *data_raw, std::string **p_string_array,
    int *p_array_length)
{
  // Check datatype version and class
  int offset = 0;
  unsigned char version_class = datatype_raw[offset++];
  if (version_class >> 4 != 1)
    throw ray_trace_exception("Error: Unexpected HDF5 datatype version.\n");
  if ((version_class & 0b00001111) != 3)
    throw ray_trace_exception("Error: Unexpected HDF5 datatype class.\n");

  // Read datatype metadata
  unsigned char class_1 = datatype_raw[offset++];
  offset += 2;

  // Read data size
  unsigned int size;
  std::memcpy(&size, datatype_raw + offset, 4);
  offset += 4;

  // Check character set
  if (class_1 >> 4 != 0)
    throw ray_trace_exception("Error: Unexpected HDF5 string encoding.\n");

  // Read and check dimensions
  unsigned long int *dims;
  int num_dims;
  read_hdf5_dataspace_dims(dataspace_raw, &dims, &num_dims);
  if (num_dims != 1)
    throw ray_trace_exception("Error: Unexpected HDF5 string encoding.\n");
  *p_array_length = static_cast<int>(dims[0]);

  // Allocate and initialize array
  *p_string_array = new std::string[*p_array_length];
  char *buffer = new char[size];
  for (int n = 0; n < *p_array_length; n++)
  {
    std::memcpy(buffer, data_raw + size * static_cast<unsigned int>(n), size);
    (*p_string_array)[n].assign(buffer, size);
  }

  // Free dimensions
  delete[] dims;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to dimensions from raw HDF5 dataspace
// Inputs:
//   dataspace_raw: raw dataspace description
// Outputs:
//   *p_dims: newly allocated and initialized array of dimensions
//   *p_num_dims: number of dimensions
// Notes:
//   Must not have permutation indices.
//   Must be run on little-endian machine

void athena_reader::read_hdf5_dataspace_dims(const unsigned char *dataspace_raw,
    unsigned long int **p_dims, int *p_num_dims)
{
  // Read and check metadata
  int offset = 0;
  if (dataspace_raw[offset++] != 1)
    throw ray_trace_exception("Error: Unexpected HDF5 dataspace version.\n");
  *p_num_dims = dataspace_raw[offset++];
  unsigned char flags = dataspace_raw[offset++];
  if (flags & 0b00000010)
    throw ray_trace_exception("Error: Unexpected HDF5 dataspace permutation indices.\n");
  offset += 5;

  // Allocate and extract dimensions
  *p_dims = new unsigned long int[*p_num_dims];
  std::memcpy(*p_dims, dataspace_raw + offset, static_cast<std::size_t>(8 * (*p_num_dims)));
  return;
}
