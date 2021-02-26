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
  AthenaReader(const InputReader *p_input_reader);
  AthenaReader(const AthenaReader &source) = delete;
  AthenaReader &operator=(const AthenaReader &source) = delete;
  ~AthenaReader();

  // Input data - general
  ModelType model_type;

  // Input data - simulation and plasma parameters
  std::string simulation_file;
  std::string simulation_kappa_name;
  PlasmaModel plasma_model;

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

  // Data
  int n_3_root;
  Array<int> levels;
  Array<int> locations;
  Array<float> x1f, x2f, x3f;
  Array<float> x1v, x2v, x3v;
  Array<float> prim, bb;

  // External function
  void Read();

  // Internal functions
  void VerifyVariables();
  void ReadHDF5Superblock();
  void ReadHDF5RootGroupSymbolTableEntry();
  void ReadHDF5RootHeap();
  void ReadHDF5RootObjectHeader();
  void ReadHDF5Tree();
  unsigned long int ReadHDF5DatasetHeaderAddress(const char *name);
  void ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
      unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw);
  static void ReadHDF5DataspaceDims(const unsigned char *dataspace_raw, unsigned long int **p_dims,
      int *p_num_dims);
  void ReadHDF5IntArray(const char *name, Array<int> &int_array);
  void ReadHDF5FloatArray(const char *name, Array<float> &float_array);
  static void SetHDF5StringArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, std::string **string_array,
      int *p_array_length);
  static void SetHDF5IntArray(const unsigned char *datatype_raw, const unsigned char *dataspace_raw,
      const unsigned char *data_raw, Array<int> &int_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array);
};

#endif
