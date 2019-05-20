// Ray Trace Athena++ reader

// C++ headers
#include <fstream>  // ifstream
#include <string>   // string

// Ray Trace headers
#include "read_athena.hpp"
#include "read_input.hpp"   // input_reader
#include "exceptions.hpp"   // ray_trace_exception

//--------------------------------------------------------------------------------------------------

// Athena++ reader constructor
// Inputs:
//   input_file: object containing input file data
athena_reader::athena_reader(const std::string data_file_)
  : data_file(data_file_) {}

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
  // Open data file
  std::ifstream data_stream(data_file);
  if (not data_stream.is_open())
    throw ray_trace_exception("Error: Could not open data file.\n");

  // Read superblock
  unsigned long int root_object_header_address, btree_address, root_name_heap_address;
  read_hdf5_superblock(data_stream, &root_object_header_address, &btree_address,
      &root_name_heap_address);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 superblock
// Inputs:
//   data_stream: input file stream
// Outputs:
//   *p_root_object_header_address: offset into file of root header
//   *p_btree_address: offset into file of root node
//   *p_root_name_heap_address: offset into file of root heap
// Notes:
//   Does not allow for superblock to be anywhere but beginning of file.
//   Must have superblock version 0.
//   Must have size of offsets 8.
//   Must have size of lengths 8.
void athena_reader::read_hdf5_superblock(std::ifstream &data_stream,
    unsigned long int *p_root_object_header_address, unsigned long int *p_btree_address,
    unsigned long int *p_root_name_heap_address)
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
  read_root_group_symbol_table_entry(data_stream, p_root_object_header_address, p_btree_address,
      p_root_name_heap_address);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root group symbol table entry
// Inputs:
//   data_stream: input file stream
// Outputs:
//   *p_root_object_header_address: offset into file of root header
//   *p_btree_address: offset into file of root node
//   *p_root_name_heap_address: offset into file of root heap
// Notes:
//   Assumes stream pointer is already set.
//   Must have size of offsets 8.
//   Must be run on little-endian machine
void athena_reader::read_root_group_symbol_table_entry(std::ifstream &data_stream,
    unsigned long int *p_root_object_header_address, unsigned long int *p_btree_address,
    unsigned long int *p_root_name_heap_address)
{
  // Skip reading link name offset
  data_stream.ignore(8);

  // Read object header address
  data_stream.read(reinterpret_cast<char *>(p_root_object_header_address), 8);

  // Check cache type
  unsigned int cache_type;
  data_stream.read(reinterpret_cast<char *>(&cache_type), 4);
  if (cache_type != 1)
    throw ray_trace_exception("Error: Unexpected HDF5 root group symbol table entry cache type.\n");

  // Read B-tree and name heap addresses
  data_stream.read(reinterpret_cast<char *>(p_btree_address), 8);
  data_stream.read(reinterpret_cast<char *>(p_root_name_heap_address), 8);
  return;
}
