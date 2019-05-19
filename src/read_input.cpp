// Ray Trace input reader

// C++ headers
#include <filesystem>  // filesystem
#include <iostream>    // cout

// Ray Trace headers
#include "read_input.hpp"

// Input reader
// Inputs:
//   input_file: C string containing name of input file
// Outputs: (none)

void read_input(const char *input_file)
{
  unsigned long int input_file_length = std::filesystem::file_size(input_file);
  std::cout << input_file_length << "\n";
  return;
}
