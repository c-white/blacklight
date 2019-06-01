// Ray Trace Athena++ reader header

#ifndef READ_ATHENA_H_
#define READ_ATHENA_H_

// C++ headers
#include <fstream>  // ifstream
#include <string>   // string

// Ray Trace headers
#include "array.hpp"       // array
#include "read_input.hpp"  // input_reader

//--------------------------------------------------------------------------------------------------

// Athena++ reader
struct athena_reader
{
  // Constructors and destructor
  athena_reader(const std::string data_file);
  athena_reader(const athena_reader &source) = delete;
  athena_reader &operator=(const athena_reader &source) = delete;
  ~athena_reader();

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
  array<int> num_variables;
  unsigned long int *children_addresses;
  int num_children = 0;

  // Data
  array<int> levels;
  array<int> locations;
  array<float> rf, thf, phf;
  array<float> r, th, ph;
  array<float> prim, bb;

  // Functions
  void read();
  void read_hdf5_superblock();
  void read_root_group_symbol_table_entry();
  void read_hdf5_root_heap();
  void read_hdf5_root_object_header();
  void read_hdf5_tree();
  void read_hdf5_int_array(const char *name, array<int> &int_array);
  void read_hdf5_float_array(const char *name, array<float> &float_array);
  unsigned long int read_hdf5_dataset_header_address(const char *name);
  void read_hdf5_data_object_header(unsigned long int data_object_header_address,
      unsigned char **p_datatype_raw, unsigned char **p_dataspace_raw, unsigned char **p_data_raw);
  static void set_hdf5_string_array(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      std::string **string_array, int *p_array_length);
  static void set_hdf5_int_array(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, array<int> &int_array);
  static void set_hdf5_float_array(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, array<float> &float_array);
  static void read_hdf5_dataspace_dims(const unsigned char *dataspace_raw,
      unsigned long int **p_dims, int *p_num_dims);
};

#endif
