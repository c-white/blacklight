// Blacklight Athena++ reader

// C++ headers
#include <cstring>   // memcpy, size_t
#include <fstream>   // ifstream
#include <ios>       // ios_base, streamoff
#include <optional>  // optional
#include <string>    // getline, string

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "read_athena.hpp"
#include "array.hpp"        // Array
#include "exceptions.hpp"   // BlacklightException

//--------------------------------------------------------------------------------------------------

// Athena++ reader constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
// Notes:
//   Opens stream for reading.
AthenaReader::AthenaReader(const InputReader *p_input_reader)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();

  // Proceed only if needed
  if (model_type == ModelType::simulation)
  {
    // Copy simulation parameters
    simulation_file = p_input_reader->simulation_file.value();

    // Open file
    data_stream = std::ifstream(simulation_file, std::ios_base::in | std::ios_base::binary);
    if (not data_stream.is_open())
      throw BlacklightException("Could not open data file.");
  }
}

//--------------------------------------------------------------------------------------------------

// Athena++ reader destructor
AthenaReader::~AthenaReader()
{
  // Free memory
  if (num_dataset_names > 0)
    delete[] dataset_names;
  if (num_variable_names > 0)
    delete[] variable_names;
  if (num_children > 0)
    delete[] children_addresses;
}

//--------------------------------------------------------------------------------------------------

// Athena++ reader read and initialize function
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Does nothing if model does not need to be read from file.
//   Initializes all member objects.
//   Implements a subset of the HDF5 standard:
//       portal.hdfgroup.org/display/HDF5/File+Format+Specification
void AthenaReader::Read()
{
  // Only proceed if needed
  if (model_type != ModelType::simulation)
    return;

  // Read basic data about file
  ReadHDF5Superblock();
  ReadHDF5RootHeap();
  ReadHDF5RootObjectHeader();
  ReadHDF5Tree();

  // Read block layout
  ReadHDF5IntArray("Levels", levels);
  ReadHDF5IntArray("LogicalLocations", locations);

  // Read coordinates
  ReadHDF5FloatArray("x1f", x1f);
  ReadHDF5FloatArray("x2f", x2f);
  ReadHDF5FloatArray("x3f", x3f);
  ReadHDF5FloatArray("x1v", x1v);
  ReadHDF5FloatArray("x2v", x2v);
  ReadHDF5FloatArray("x3v", x3v);

  // Read cell data
  ReadHDF5FloatArray("prim", prim);
  ReadHDF5FloatArray("B", bb);

  // Close data file
  data_stream.close();
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
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, root_grid_size);
        n_3_root = root_grid_size(2);
      }
      else if (name == "DatasetNames")
      {
        dataset_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset, &dataset_names,
            &num_dataset_names);
      }
      else if (name == "VariableNames")
      {
        variable_names_found = true;
        SetHDF5StringArray(datatype_raw, dataspace_raw, message_data + offset, &variable_names,
            &num_variable_names);
      }
      else if (name == "NumVariables")
      {
        num_variables_found = true;
        SetHDF5IntArray(datatype_raw, dataspace_raw, message_data + offset, num_variables);
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

// Function to check that needed variables located as expected
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes metadata set.
void AthenaReader::VerifyVariables()
{
  // Check that primitives and magnetic field are present in expected order
  if (num_dataset_names < 2)
    throw BlacklightException("Insufficient datasets in data file.");
  if (not dataset_names[ind_prim].compare("prim"))
    throw BlacklightException("Primitives not found in data file.");
  if (not dataset_names[ind_bb].compare("B"))
    throw BlacklightException("Magnetic field not found in data file.");

  // Check that primitives and magnetic field have expected sizes
  if (num_variables(ind_prim) < 5)
    throw BlacklightException("Primitives from data file do not have 5 variables.");
  if (num_variables(ind_bb) < 3)
    throw BlacklightException("Magnetic field from data file does not have 3 variables.");

  // Check that variables are in expected locations
  if (not variable_names[ind_rho].compare("rho"))
    throw BlacklightException("Density not found.");
  if (not variable_names[ind_pgas].compare("press"))
    throw BlacklightException("Gas pressure not found.");
  if (not variable_names[ind_uu1].compare("vel1"))
    throw BlacklightException("x1-velocity not found.");
  if (not variable_names[ind_uu2].compare("vel2"))
    throw BlacklightException("x2-velocity not found.");
  if (not variable_names[ind_uu3].compare("vel3"))
    throw BlacklightException("x3-velocity not found.");
  if (not variable_names[5+ind_bb1].compare("Bcc1"))
    throw BlacklightException("x1-field not found.");
  if (not variable_names[5+ind_bb2].compare("Bcc2"))
    throw BlacklightException("x2-field not found.");
  if (not variable_names[5+ind_bb3].compare("Bcc3"))
    throw BlacklightException("x3-field not found.");
  return;
}

//--------------------------------------------------------------------------------------------------

// Function to locate children of root node in HDF5 tree
// Inputs: (none)
// Outputs: (none)
//   returned value: address within file of header
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
  num_children = num_entries;

  // Skip addresses of siblings
  data_stream.ignore(16);

  // Read addresses of children
  children_addresses = new unsigned long int[num_children];
  for (int n = 0; n < num_children; n++)
  {
    data_stream.ignore(8);
    data_stream.read(reinterpret_cast<char *>(children_addresses + n), 8);
  }
  return;
}

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
//   int_array: array allocated and set (indirectly)
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

// Function to locate HDF5 header address for dataset with given name
// Inputs:
//   name: name of dataset
// Outputs:
//   returned value: address within file of header
// Notes:
//   Assumes children_addresses set.
//   Changes stream pointer.
//   Must have symbol table entry cache type 0.
//   Must have size of offsets 8.
//   Must be run on little-endian machine.
unsigned long int AthenaReader::ReadHDF5DatasetHeaderAddress(const char *name)
{
  // Go through children
  for (int n = 0; n < num_children; n++)
  {
    // Check symbol table node signature and version
    data_stream.seekg(static_cast<std::streamoff>(children_addresses[n]));
    const unsigned char expected_signature[] = {'S', 'N', 'O', 'D'};
    for (int m = 0; m < 4; m++)
      if (data_stream.get() != expected_signature[m])
        throw BlacklightException("Unexpected HDF5 symbol table node signature.");
    if (data_stream.get() != 1)
      throw BlacklightException("Unexpected HDF5 symbol table node version.");
    data_stream.ignore(1);

    // Read number of symbols
    unsigned short int num_symbols;
    data_stream.read(reinterpret_cast<char *>(&num_symbols), 2);

    // Go through symbols
    for (int m = 0; m < num_symbols; m++)
    {
      // Read addresses
      unsigned long int link_name_offset, object_header_address;
      data_stream.read(reinterpret_cast<char *>(&link_name_offset), 8);
      data_stream.read(reinterpret_cast<char *>(&object_header_address), 8);

      // Check cache type
      unsigned int cache_type;
      data_stream.read(reinterpret_cast<char *>(&cache_type), 4);
      if (cache_type != 0)
        throw BlacklightException("Unexpected HDF5 symbol table entry cache type.");

      // Skip remaining entry
      data_stream.ignore(20);

      // Compare name
      std::streamoff position = data_stream.tellg();
      data_stream.seekg(static_cast<std::streamoff>(root_data_segment_address + link_name_offset));
      std::string dataset_name;
      std::getline(data_stream, dataset_name, '\0');
      if (dataset_name == name)
        return object_header_address;
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
void AthenaReader::ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
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

// Function to initialize string array dataset from HDF5
// Inputs:
//   datatype_raw: raw datatype description
//   dataspace_raw: raw dataspace description
//   data_raw: raw data
// Outputs:
//   *p_string_array: array allocated and set
//   *p_array_length: number of allocated members
// Notes:
//   Must have datatype version 1.
//   Must be a 1D array.
void AthenaReader::SetHDF5StringArray(const unsigned char *datatype_raw,
    const unsigned char *dataspace_raw, const unsigned char *data_raw, std::string **p_string_array,
    int *p_array_length)
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
  *p_array_length = static_cast<int>(dims[0]);

  // Allocate and initialize array
  *p_string_array = new std::string[*p_array_length];
  char *buffer = new char[size];
  for (int n = 0; n < *p_array_length; n++)
  {
    std::memcpy(buffer, data_raw + size * static_cast<unsigned int>(n), size);
    (*p_string_array)[n].assign(buffer, size);
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
//   int_array: array allocated and set
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
    int_array.Allocate(static_cast<int>(dims[0]));
  else if (num_dims == 2)
    int_array.Allocate(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
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
//   float_array: array allocated and set
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
  if (num_dims == 1)
    float_array.Allocate(static_cast<int>(dims[0]));
  else if (num_dims == 2)
    float_array.Allocate(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
  else if (num_dims == 5)
    float_array.Allocate(static_cast<int>(dims[0]), static_cast<int>(dims[1]),
        static_cast<int>(dims[2]), static_cast<int>(dims[3]), static_cast<int>(dims[4]));
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
void AthenaReader::ReadHDF5DataspaceDims(const unsigned char *dataspace_raw,
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
  *p_dims = new unsigned long int[*p_num_dims];
  std::memcpy(*p_dims, dataspace_raw + offset, static_cast<std::size_t>(8 * (*p_num_dims)));
  return;
}
