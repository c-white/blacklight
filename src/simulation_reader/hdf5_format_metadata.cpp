// Blacklight simulation reader - HDF5 interface for reading metadata

// C++ headers
#include <cstring>  // memcpy, size_t
#include <fstream>  // ifstream
#include <ios>      // streamoff
#include <string>   // getline, string

// Blacklight headers
#include "simulation_reader.hpp"
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function to locate HDF5 header address for dataset with given name
// Inputs:
//   name: null-terminated name of dataset, possibly including intermediate slashes, excluding
//       leading slash
//   btree_address: address of B-tree node containing dataset
//   data_segment_address: address of data segment in local name heap of B-tree
// Outputs:
//   returned value: address within file of header
// Notes:
//   Changes stream pointer.
//   Must have B-tree version 1.
//   Must have symbol table entry cache type 0 or 1.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
unsigned long int SimulationReader::ReadHDF5DatasetHeaderAddress(const char *name,
    unsigned long int btree_address, unsigned long int data_segment_address)
{
  // Parse name
  std::string name_str(name);
  std::string::size_type slash_pos = name_str.find("/");

  // Check tree signature
  data_stream.seekg(static_cast<std::streamoff>(btree_address));
  const unsigned char expected_tree_signature[] = {'T', 'R', 'E', 'E'};
  for (int n = 0; n < 4; n++)
    if (data_stream.get() != expected_tree_signature[n])
      throw BlacklightException("Unexpected HDF5 B-tree signature.");

  // Check node type
  if (data_stream.get() != 0)
    throw BlacklightException("Unexpected HDF5 node type.");

  // Skip node level
  data_stream.ignore(1);

  // Read number of children
  unsigned short int num_children;
  data_stream.read(reinterpret_cast<char *>(&num_children), 2);

  // Skip addresses of siblings
  data_stream.ignore(16);

  // Go through children
  std::streamoff list_begin = data_stream.tellg();
  for (int n_child = 0; n_child < num_children; n_child++)
  {
    // Get child address
    data_stream.seekg(list_begin + 16 * n_child + 8);
    unsigned long int child_address = 0;
    data_stream.read(reinterpret_cast<char *>(&child_address), 8);

    // Check symbol table node signature and version
    data_stream.seekg(static_cast<std::streamoff>(child_address));
    const unsigned char expected_symbol_table_signature[] = {'S', 'N', 'O', 'D'};
    for (int m = 0; m < 4; m++)
      if (data_stream.get() != expected_symbol_table_signature[m])
        throw BlacklightException("Unexpected HDF5 symbol table node signature.");
    if (data_stream.get() != 1)
      throw BlacklightException("Unexpected HDF5 symbol table node version.");
    data_stream.ignore(1);

    // Read number of symbols
    unsigned short int num_symbols;
    data_stream.read(reinterpret_cast<char *>(&num_symbols), 2);

    // Go through symbols
    for (int n_symbol = 0; n_symbol < num_symbols; n_symbol++)
    {
      // Read addresses
      unsigned long int link_name_offset, object_header_address;
      data_stream.read(reinterpret_cast<char *>(&link_name_offset), 8);
      data_stream.read(reinterpret_cast<char *>(&object_header_address), 8);

      // Check cache type
      unsigned int cache_type;
      data_stream.read(reinterpret_cast<char *>(&cache_type), 4);
      data_stream.ignore(4);
      if (cache_type == 0 or cache_type == 1)
        data_stream.ignore(16);
      else
        throw BlacklightException("Unexpected HDF5 symbol table entry cache type.");

      // Search for match to dataset name
      std::streamoff position = data_stream.tellg();
      data_stream.seekg(static_cast<std::streamoff>(data_segment_address + link_name_offset));
      std::string local_name;
      std::getline(data_stream, local_name, '\0');
      if (slash_pos == std::string::npos and local_name == name_str)
        return object_header_address;

      // Search for match to group name
      if (slash_pos != std::string::npos and name_str.compare(0, slash_pos, local_name) == 0)
      {
        // Go to data object header
        data_stream.seekg(static_cast<std::streamoff>(object_header_address));

        // Check data object header version
        if (data_stream.get() != 1)
          throw BlacklightException("Unexpected HDF5 object header version.");
        data_stream.ignore(1);

        // Get number of messages
        unsigned short int num_messages;
        data_stream.read(reinterpret_cast<char *>(&num_messages), 2);

        // Skip to list of messages
        data_stream.ignore(12);

        // Go through messages
        for (int n_message = 0; n_message < num_messages; n_message++)
        {
          unsigned short int message_type;
          data_stream.read(reinterpret_cast<char *>(&message_type), 2);
          if (message_type != 17)
            continue;
          unsigned short int message_size;
          data_stream.read(reinterpret_cast<char *>(&message_size), 2);
          if (message_size != 16)
            throw BlacklightException("Unexpected HDF5 header message size.");
          unsigned char message_flags;
          data_stream.read(reinterpret_cast<char *>(&message_flags), 1);
          if (message_flags & 0b00000010)
            throw BlacklightException("Unexpected HDF5 header message flag.");
          data_stream.ignore(3);
          unsigned long int btree_node_address, heap_node_address;
          data_stream.read(reinterpret_cast<char *>(&btree_node_address), 8);
          data_stream.read(reinterpret_cast<char *>(&heap_node_address), 8);
          unsigned long int data_segment_node_address = ReadHDF5Heap(heap_node_address);
          return ReadHDF5DatasetHeaderAddress(name + slash_pos + 1, btree_node_address,
              data_segment_node_address);
        }
      }

      // Prepare for next symbol
      data_stream.seekg(position);
    }
  }

  // Report failure to find named dataset
  throw BlacklightException("Could not find HDF5 dataset in file.");
  return 0;
}

//--------------------------------------------------------------------------------------------------

// Function to read HDF5 data object header
// Inputs:
//   data_object_header_address: offset where header is located
// Outputs:
//   *p_datatype_raw: raw datatype description
//   *p_dataspace_raw: raw dataspace description
//   *p_data_raw: raw data
// Notes:
//   Changes stream pointer.
//   Must have object header version 1.
//   Must not have shared header messages.
//   Must have data layout message version 3.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
void SimulationReader::ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
    unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw)
{
  // Check object header version
  data_stream.seekg(static_cast<std::streamoff>(data_object_header_address));
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

  // Prepare containers for message information
  bool datatype_found = false;
  bool dataspace_found = false;
  bool data_layout_found = false;
  unsigned long int data_address = 0;
  unsigned long int data_size = 0;

  // Go through messages
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

    // Inspect any datatype messages
    }
    else if (message_type == 3)
    {
      if (datatype_found)
        throw BlacklightException("Too many HDF5 datatypes for dataset.");
      datatype_found = true;
      *p_datatype_raw = new unsigned char[message_size];
      std::memcpy(*p_datatype_raw, message_data, message_size);
    }

    // Inspect any dataspace messages
    else if (message_type == 1)
    {
      if (dataspace_found)
        throw BlacklightException("Too many HDF5 dataspaces for dataset.");
      dataspace_found = true;
      *p_dataspace_raw = new unsigned char[message_size];
      std::memcpy(*p_dataspace_raw, message_data, message_size);
    }

    // Inspect any data layout messages
    else if (message_type == 8)
    {

      // Note if data layout has already been found
      if (data_layout_found)
        throw BlacklightException("Too many HDF5 data layouts for dataset.");
      data_layout_found = true;

      // Check layout version and class
      int offset = 0;
      if (message_data[offset++] != 3)
        throw BlacklightException("Unexpected HDF5 data layout message version.");
      if (message_data[offset++] != 1)
        throw BlacklightException("Unexpected HDF5 data layout class.");

      // Read data layout
      std::memcpy(&data_address, message_data + offset, 8);
      offset += 8;
      std::memcpy(&data_size, message_data + offset, 8);
      offset += 8;
    }

    // Free raw buffer
    delete[] message_data;

    // Break when required information found
    if (datatype_found and dataspace_found and data_layout_found)
      break;
  }

  // Check that appropriate messages were found
  if (not (datatype_found and dataspace_found and data_layout_found))
    throw BlacklightException("Could not find needed dataset properties.");

  // Read raw data
  *p_data_raw = new unsigned char[data_size];
  data_stream.seekg(static_cast<std::streamoff>(data_address));
  data_stream.read(reinterpret_cast<char *>(*p_data_raw), static_cast<std::streamoff>(data_size));
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to read dimensions from raw HDF5 dataspace
// Inputs:
//   dataspace_raw: raw dataspace description
// Outputs:
//   *p_dims: newly allocated and initialized array of dimensions
//   *p_num_dims: number of dimensions
// Notes:
//   Must not have permutation indices.
//   Must be run on little-endian machine.
void SimulationReader::ReadHDF5DataspaceDims(const unsigned char *dataspace_raw,
    unsigned long int **p_dims, int *p_num_dims)
{
  // Read and check metadata
  int offset = 0;
  if (dataspace_raw[offset++] != 1)
    throw BlacklightException("Unexpected HDF5 dataspace version.");
  *p_num_dims = dataspace_raw[offset++];
  unsigned char flags = dataspace_raw[offset++];
  if (flags & 0b00000010)
    throw BlacklightException("Unexpected HDF5 dataspace permutation indices.");
  offset += 5;

  // Allocate and extract dimensions
  if (*p_num_dims > 0)
  {
    *p_dims = new unsigned long int[*p_num_dims];
    std::memcpy(*p_dims, dataspace_raw + offset, static_cast<std::size_t>(8 * (*p_num_dims)));
  }
  else
    *p_dims = nullptr;
  return;
}
