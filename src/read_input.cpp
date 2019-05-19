// Ray Trace input reader

// C++ headers
#include <filesystem>  // filesystem
#include <iostream>    // cout

// Ray Trace headers
#include "read_input.hpp"

// Input reader constructor
// Inputs:
//   input_file: C string containing name of input file

input_reader::input_reader(const char *input_file_)
  : input_file(input_file_)
{
  unsigned long int input_file_length = std::filesystem::file_size(input_file);
  std::cout << input_file_length << "\n";
  return;
}

// Input reader destructor

input_reader::~input_reader() {}
