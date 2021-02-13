// Blacklight main file

// C++ headers
#include <iostream>  // cout
#include <optional>  // bad_optional_access, optional
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
  int num_runs;
  InputReader *p_input_reader;
  try
  {
    p_input_reader = new InputReader(input_file);
    num_runs = p_input_reader->Read();
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
  try
  {
    omp_set_num_threads(p_input_reader->num_threads.value());
  }
  catch (const std::bad_optional_access &exception)
  {
    std::cout << "Error: num_threads not specified in input file.\n";
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not set number of threads.\n";
    return 1;
  }

  // Go through runs
  for (int n = 0; n < num_runs; n++)
  {
    // Adjust file names
    p_input_reader->AdjustFileNames(n);

    // Read data file
    AthenaReader *p_athena_reader;
    try
    {
      p_athena_reader = new AthenaReader(p_input_reader);
      p_athena_reader->Read();
    }
    catch (const BlacklightException &exception)
    {
      std::cout << exception.what();
      return 1;
    }
    catch (const std::bad_optional_access &exception)
    {
      std::cout << "Error: AthenaReader unable to find all needed values in input file.\n";
      return 1;
    }
    catch (...)
    {
      std::cout << "Error: Could not read data file.\n";
      return 1;
    }

    // Process data
    RayTracer *p_ray_tracer;
    try
    {
      p_ray_tracer = new RayTracer(p_input_reader, p_athena_reader);
      p_ray_tracer->MakeImage();
    }
    catch (const BlacklightException &exception)
    {
      std::cout << exception.what();
      return 1;
    }
    catch (const std::bad_optional_access &exception)
    {
      std::cout << "Error: RayTracer unable to find all needed values in input file.\n";
      return 1;
    }
    catch (...)
    {
      std::cout << "Error: Could not process data.\n";
      return 1;
    }

    // Write output file
    OutputWriter *p_output_writer;
    try
    {
      p_output_writer = new OutputWriter(p_input_reader, p_ray_tracer);
      p_output_writer->Write();
    }
    catch (const BlacklightException &exception)
    {
      std::cout << exception.what();
      return 1;
    }
    catch (const std::bad_optional_access &exception)
    {
      std::cout << "Error: OutputWriter unable to find all needed values in input file.\n";
      return 1;
    }
    catch (...)
    {
      std::cout << "Error: Could not write output file.\n";
      return 1;
    }

    // Free memory
    delete p_output_writer;
    delete p_ray_tracer;
    delete p_athena_reader;
  }

  // Free memory
  delete p_input_reader;

  // End program
  return 0;
}
