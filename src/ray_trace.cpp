// Ray Trace main file

// C++ headers
#include <filesystem>  // filesystem
#include <iostream>    // cout

// Ray Trace headers
#include "ray_trace.hpp"
#include "read_input.hpp"  // read_input

// Main function
int main(int argc, char *argv[])
{
  // Parse command line inputs
  if (argc != 2)
  {
    std::cout << "Error: Must give a single input file.\n";
    return 1;
  }
  const char *input_file = argv[1];

  // Read input file
  read_input(input_file);

  // End program
  return 0;
}
