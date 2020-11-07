// Blacklight main file

// C++ headers
#include <iostream>  // cout
#include <string>    // string

// Library headers
#include <omp.h>  // omp_set_num_threads

// Blacklight headers
#include "blacklight.hpp"
#include "exceptions.hpp"    // BlacklightException
#include "ray_tracer.hpp"    // RayTracer
#include "read_athena.hpp"   // AthenaReader
#include "read_input.hpp"    // InputReader
#include "write_output.hpp"  // OutputReader

//--------------------------------------------------------------------------------------------------

// Main function
// Inputs:
//   argc: number of command-line arguments
//   argv: array of command-line argument strings
// Outputs:
//   returned value:
//     0: successful execution
//     1: error
int main(int argc, char *argv[])
{
  // Parse command line inputs
  if (argc != 2)
  {
    std::cout << "Error: Must give a single input file.\n";
    return 1;
  }
  const std::string input_file(argv[1]);

  // Read input file
  InputReader input_reader(input_file);
  try
  {
    input_reader.Read();
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not read input file.\n";
    return 1;
  }

  // Set number of threads to use
  omp_set_num_threads(input_reader.num_threads);

  // Read data file
  AthenaReader athena_reader(input_reader);
  try
  {
    athena_reader.Read();
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not read data file.\n";
    return 1;
  }

  // Process data
  RayTracer ray_tracer(input_reader, athena_reader);
  try
  {
    ray_tracer.MakeImage();
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not process data.\n";
    return 1;
  }

  // Write output file
  OutputWriter output_writer(input_reader, ray_tracer);
  try
  {
    output_writer.Write();
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not write output file.\n";
    return 1;
  }

  // End program
  return 0;
}
