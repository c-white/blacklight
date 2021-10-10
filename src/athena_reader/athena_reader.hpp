// Blacklight Athena++ reader header

#ifndef ATHENA_READER_H_
#define ATHENA_READER_H_

// C++ headers
#include <fstream>  // ifstream
#include <string>   // string

// Blacklight headers
#include "../blacklight.hpp"                 // enums
#include "../input_reader/input_reader.hpp"  // InputReader
#include "../utils/array.hpp"                // Array

//--------------------------------------------------------------------------------------------------

// Athena++ reader
struct AthenaReader
{
  // Constructors and destructor
  AthenaReader(const InputReader *p_input_reader_);
  AthenaReader(const AthenaReader &source) = delete;
  AthenaReader &operator=(const AthenaReader &source) = delete;
  ~AthenaReader();

  // Pointers to other objects
  const InputReader *p_input_reader;

  // Input data - general
  ModelType model_type;

  // Input data - simulation parameters
  std::string simulation_file;
  bool simulation_multiple;
  int simulation_start;
  int simulation_end;
  std::string simulation_kappa_name;

  // Input data - plasma parameters
  PlasmaModel plasma_model;

  // Input data - slow light parameters
  bool slow_light_on;
  int slow_chunk_size;
  double slow_t_start;
  double slow_dt;

  // Flags for tracking function calls
  bool first_time = true;
  bool first_time_root_object_header = true;
  bool first_time_tree = true;

  // Metadata
  std::ifstream data_stream;
  unsigned long int root_object_header_address;
  unsigned long int btree_address;
  unsigned long int root_name_heap_address;
  unsigned long int root_data_segment_address;
  std::string *dataset_names;
  int num_dataset_names = 0;
  std::string *variable_names;
  int num_variable_names = 0;
  Array<int> num_variables;
  unsigned long int *children_addresses;
  int num_children = 0;
  int ind_rho, ind_pgas, ind_kappa;
  int ind_uu1, ind_uu2, ind_uu3;
  int ind_bb1, ind_bb2, ind_bb3;
  int num_arrays;
  int latest_file_number;
  const float extrapolation_tolerance = 1.0f;

  // Data
  int n_3_root;
  Array<int> levels;
  Array<int> locations;
  Array<float> x1f, x2f, x3f;
  Array<float> x1v, x2v, x3v;
  float *time;
  Array<float> *prim;
  Array<float> *bb;

  // External function
  double Read(int snapshot);

  // Internal functions - athena_reader.cpp
  std::string FormatFilename(int file_number);
  void VerifyVariables();

  // Internal functions - hdf5_format_structure.cpp
  void ReadHDF5Superblock();
  void ReadHDF5RootGroupSymbolTableEntry();
  void ReadHDF5RootHeap();
  void ReadHDF5RootObjectHeader();
  void ReadHDF5Tree();
  void ReadHDF5FloatAttribute(const char *attribute_name, float *p_val);

  // Internal functions - hdf5_format_metadata.cpp
  unsigned long int ReadHDF5DatasetHeaderAddress(const char *name);
  void ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
      unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw);
  static void ReadHDF5DataspaceDims(const unsigned char *dataspace_raw, unsigned long int **p_dims,
      int *p_num_dims);

  // Internal functions - hdf5_format_arrays.cpp
  void ReadHDF5IntArray(const char *name, Array<int> &int_array);
  void ReadHDF5FloatArray(const char *name, Array<float> &float_array);
  static void SetHDF5StringArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, bool allocate,
      std::string **string_array, int *p_array_length);
  static void SetHDF5IntArray(const unsigned char *datatype_raw, const unsigned char *dataspace_raw,
      const unsigned char *data_raw, Array<int> &int_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array);
};

#endif
