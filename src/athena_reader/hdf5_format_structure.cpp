// Blacklight Athena++ reader - HDF5 general structure interface

// C++ headers
#include <cstring>  // memcpy
#include <fstream>  // ifstream
#include <ios>      // streamoff
#include <string>   // string

// Blacklight headers
#include "athena_reader.hpp"
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

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
void AthenaReader::ReadHDF5Superblock()
{
  // Check format signature
  data_stream.seekg(0);
  const unsigned char expected_signature[] = {0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a};
  for (int n = 0; n < 8; n++)
    if (data_stream.get() != expected_signature[n])
      throw BlacklightException("Unexpected HDF5 format signature.");

  // Check superblock version
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 superblock version.");

  // Check other version numbers
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 file free space storage version.");
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 root group symbol table entry version.");
  data_stream.ignore(1);
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 shared header message format version.");

  // Check sizes
  if (data_stream.get() != 8)
    throw BlacklightException("Unexpected HDF5 size of offsets.");
  if (data_stream.get() != 8)
    throw BlacklightException("Unexpected HDF5 size of lengths.");
  data_stream.ignore(1);

  // Skip checking tree parameters and consistency flags
  data_stream.ignore(2 * 2 + 4);

  // Skip checking addresses
  data_stream.ignore(4 * 8);

  // Read root group symbol table entry
  ReadHDF5RootGroupSymbolTableEntry();
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
//   Must be run on little-endian machine.
void AthenaReader::ReadHDF5RootGroupSymbolTableEntry()
{
  // Skip reading link name offset
  data_stream.ignore(8);

  // Read object header address
  data_stream.read(reinterpret_cast<char *>(&root_object_header_address), 8);

  // Check cache type
  unsigned int cache_type;
  data_stream.read(reinterpret_cast<char *>(&cache_type), 4);
  if (cache_type != 1)
    throw BlacklightException("Unexpected HDF5 root group symbol table entry cache type.");
  data_stream.ignore(4);

  // Read B-tree and name heap addresses
  data_stream.read(reinterpret_cast<char *>(&btree_address), 8);
  data_stream.read(reinterpret_cast<char *>(&root_name_heap_address), 8);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root local heap
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets root_data_segment_address.
//   Assumes root_name_heap_address set.
//   Changes stream pointer.
//   Must have size of offsets 8.
void AthenaReader::ReadHDF5RootHeap()
{
  // Check local heap signature and version
  data_stream.seekg(static_cast<std::streamoff>(root_name_heap_address));
  const unsigned char expected_signature[] = {'H', 'E', 'A', 'P'};
  for (int n = 0; n < 4; n++)
    if (data_stream.get() != expected_signature[n])
      throw BlacklightException("Unexpected HDF5 heap signature.");
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 heap version.");
  data_stream.ignore(3);

  // Skip data segement size and offset to head of free list
  data_stream.ignore(16);

  // Read address of data segment
  data_stream.read(reinterpret_cast<char *>(&root_data_segment_address), 8);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 root object header
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Allocates and sets dataset_names, variable_names, nums_variables.
//   Sets num_dataset_names and num_variable_names.
//   Assumes root_object_header_address set.
//   Changes stream pointer.
//   Must have object header version 1.
//   Must not have shared header messages.
//   Must have attribute message version 1.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void AthenaReader::ReadHDF5RootObjectHeader()
{
  // Check object header version
  data_stream.seekg(static_cast<std::streamoff>(root_object_header_address));
  if (data_stream.get() != 1)
    throw BlacklightException("Unexpected HDF5 object header version.");
  data_stream.ignore(1);

  // Read number of header messages
  unsigned short int num_messages;
  data_stream.read(reinterpret_cast<char *>(&num_messages), 2);

  // Skip reading object reference count and object header size
  data_stream.ignore(8);

  // Align to 8 bytes within header (location of padding not documented)
  data_stream.ignore(4);

  // Go through messages
  bool root_grid_size_found = false;
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
      throw BlacklightException("Unexpected HDF5 header message flag.");

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
    }

    // Inspect any attribute messages
    else if (message_type == 12)
    {

      // Check attribute message version
      int offset = 0;
      if (message_data[offset] != 1)
        throw BlacklightException("Unexpected HDF5 attribute message version.");
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
      if (name == "RootGridSize")
      {
        root_grid_size_found = true;
        Array<int> root_grid_size;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, true, root_grid_size);
        n_3_root = root_grid_size(2);
      }
      else if (name == "DatasetNames")
      {
        dataset_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset, first_time,
            &dataset_names, &num_dataset_names);
      }
      else if (name == "VariableNames")
      {
        variable_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset, first_time,
            &variable_names, &num_variable_names);
      }
      else if (name == "NumVariables")
      {
        num_variables_found = true;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, first_time,
            num_variables);
      }

      // Free raw buffers
      delete[] name_raw;
      delete[] datatype_raw;
      delete[] dataspace_raw;
    }

    // Free raw buffer
    delete[] message_data;

    // Break when required information found
    if (root_grid_size_found and dataset_names_found and variable_names_found
        and num_variables_found)
      break;
  }

  // Check that appropriate messages were found
  if (not (root_grid_size_found and dataset_names_found and variable_names_found
      and num_variables_found))
    throw BlacklightException("Could not find needed file-level attributes.");
  if (num_variables.n1 != num_dataset_names)
    throw BlacklightException("DatasetNames and NumVariables file-level attribute mismatch.");
  VerifyVariables();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to locate children of root node in HDF5 tree
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Allocates and sets children_addresses.
//   Sets num_children.
//   Assumes btree_address set.
//   Changes stream pointer.
//   Must have B-tree version 1.
//   Must have only 1 level of children.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void AthenaReader::ReadHDF5Tree()
{
  // Check signature
  data_stream.seekg(static_cast<std::streamoff>(btree_address));
  const unsigned char expected_signature[] = {'T', 'R', 'E', 'E'};
  for (int n = 0; n < 4; n++)
    if (data_stream.get() != expected_signature[n])
      throw BlacklightException("Unexpected HDF5 B-tree signature.");

  // Check node type and level
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 node type.");
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 node level.");

  // Read number of children
  unsigned short int num_entries;
  data_stream.read(reinterpret_cast<char *>(&num_entries), 2);
  if (first_time)
    num_children = num_entries;
  else if (num_children != num_entries)
    throw BlacklightException("File layout mismatch upon subsequent read.");

  // Skip addresses of siblings
  data_stream.ignore(16);

  // Read addresses of children
  if (first_time)
    children_addresses = new unsigned long int[num_children];
  for (int n = 0; n < num_children; n++)
  {
    data_stream.ignore(8);
    data_stream.read(reinterpret_cast<char *>(children_addresses + n), 8);
  }
  return;
}
