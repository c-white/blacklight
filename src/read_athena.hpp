// Ray Trace Athena++ reader header

#ifndef READ_ATHENA_H_
#define READ_ATHENA_H_

// C++ headers
#include <string>  // string

// Ray Trace headers
#include "read_input.hpp"  // input_reader

// Input reader
struct athena_reader
{
  // Constructor and destructor
  athena_reader(const std::string data_file_);
  ~athena_reader();

  // Data
  const std::string data_file;

  // Functions
  void read();
  void read_hdf5_superblock(std::ifstream &data_stream, unsigned long int *p_root_object_header_address,
      unsigned long int *p_btree_address, unsigned long int *p_root_name_heap_address);
  void read_root_group_symbol_table_entry(std::ifstream &data_stream, unsigned long int *p_root_object_header_address,
      unsigned long int *p_btree_address, unsigned long int *p_root_name_heap_address);
};

#endif
