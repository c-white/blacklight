// Ray Trace Athena++ reader header

#ifndef READ_ATHENA_H_
#define READ_ATHENA_H_

// Ray Trace headers
#include "read_input.hpp"  // input_reader

// Input reader
struct athena_reader
{
  // Constructor and destructor
  athena_reader(const input_reader &inputs);
  ~athena_reader();

  // Data
  double r_min;
};

#endif
