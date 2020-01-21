// Ray Trace Athena++ reader header

#ifndef READ_ATHENA_H_
#define READ_ATHENA_H_

// C++ headers
#include <fstream>  // ifstream
#include <string>   // string

// Ray Trace headers
#include "array.hpp"  // Array

//--------------------------------------------------------------------------------------------------

// Athena++ reader
struct AthenaReader
{
  // Constructors and destructor
  AthenaReader(const std::string data_file);
  AthenaReader(const AthenaReader &source) = delete;
  AthenaReader &operator=(const AthenaReader &source) = delete;
  ~AthenaReader();

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
  const int ind_prim = 0;
  const int ind_rho = 0;
  const int ind_pgas = 1;
  const int ind_uu1 = 2;
  const int ind_uu2 = 3;
  const int ind_uu3 = 4;
  const int ind_bb = 1;
  const int ind_bb1 = 0;
  const int ind_bb2 = 1;
  const int ind_bb3 = 2;

  // Data
  float x1_min, x1_max, x2_min, x2_max, x3_min, x3_max;
  Array<int> levels;
  Array<int> locations;
  Array<float> x1f, x2f, x3f;
  Array<float> x1, x2, x3;
  Array<float> prim, bb;

  // Functions
  void Read();
  void ReadHDF5Superblock();
  void ReadRootGroupSymbolTableEntry();
  void ReadHDF5RootHeap();
  void ReadHDF5RootObjectHeader();
  void VerifyVariables();
  void ReadHDF5Tree();
  void ReadHDF5IntArray(const char *name, Array<int> &int_array);
  void ReadHDF5FloatArray(const char *name, Array<float> &float_array);
  unsigned long int ReadHDF5DatasetHeaderAddress(const char *name);
  void ReadHDF5DataObjectHeader(unsigned long int data_object_header_address,
      unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw);
  static void SetHDF5StringArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, std::string **string_array,
      int *p_array_length);
  static void SetHDF5IntArray(const unsigned char *datatype_raw, const unsigned char *dataspace_raw,
      const unsigned char *data_raw, Array<int> &int_array);
  static void SetHDF5FloatArray(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, Array<float> &float_array);
  static void ReadHDF5DataspaceDims(const unsigned char *dataspace_raw, unsigned long int **p_dims,
      int *p_num_dims);
};

#endif
