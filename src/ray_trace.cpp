// Ray Trace main file

// C++ headers
#include <iostream>   // cout
#include <string>     // string

// Ray Trace headers
#include "ray_trace.hpp"
#include "read_athena.hpp"  // athena_reader
#include "read_input.hpp"   // input_reader
#include "exceptions.hpp"   // ray_trace_exception

// Main function
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
  input_reader *p_inputs;
  try
  {
    p_inputs = new input_reader(input_file);
  } catch (const ray_trace_exception &exception) {
    std::cout << exception.what();
    return 1;
  } catch (...) {
    std::cout << "Error: Could not read input file.\n";
    return 1;
  }

  // Read data file
  athena_reader *p_raw_data;
  try
  {
    p_raw_data = new athena_reader(*p_inputs);
  } catch (const ray_trace_exception &exception) {
    std::cout << exception.what();
    return 1;
  } catch (...) {
    std::cout << "Error: Could not read data file.\n";
    return 1;
  }

  // End program
  return 0;
}
