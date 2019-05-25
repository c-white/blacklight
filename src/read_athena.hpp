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

  // Data
  std::ifstream data_stream;
  unsigned long int root_object_header_address;
  unsigned long int btree_address;
  unsigned long int root_name_heap_address;
  std::string *dataset_names;
  int num_dataset_names;
  std::string *variable_names;
  int num_variable_names;
  array<int> num_variables;

  // Functions
  void read();
  void read_hdf5_superblock();
  void read_root_group_symbol_table_entry();
  void read_hdf5_root_object_header();
  static void set_hdf5_string_array(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw,
      std::string **string_array, int *p_array_length);
  static void set_hdf5_int_array(const unsigned char *datatype_raw,
      const unsigned char *dataspace_raw, const unsigned char *data_raw, array<int> &int_array);
  static void read_hdf5_dataspace_dims(const unsigned char *dataspace_raw,
      unsigned long int **p_dims, int *p_num_dims);
};

#endif
