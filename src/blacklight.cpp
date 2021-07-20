// Blacklight main file

// C++ headers
#include <iostream>  // cout
#include <optional>  // bad_optional_access, optional
#include <string>    // string

// Library headers
#include <omp.h>  // omp_set_num_threads

// Blacklight headers
#include "blacklight.hpp"
#include "athena_reader/athena_reader.hpp"                // AthenaReader
#include "geodesic_integrator/geodesic_integrator.hpp"    // GeodesicIntegrator
#include "input_reader/input_reader.hpp"                  // InputReader
#include "output_writer/output_writer.hpp"                // OutputWriter
#include "radiation_integrator/radiation_integrator.hpp"  // RadiationIntegrator
#include "utils/exceptions.hpp"                           // BlacklightException

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
  // Parse command-line inputs
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

  // Define camera and integrate geodesics
  GeodesicIntegrator *p_geodesic_integrator;
  try
  {
    p_geodesic_integrator = new GeodesicIntegrator(p_input_reader);
    p_geodesic_integrator->Integrate();
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (const std::bad_optional_access &exception)
  {
    std::cout << "Error: GeodesicIntegrator unable to find all needed values in input file.\n";
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not integrate geodesics.\n";
    return 1;
  }

  // Set up Athena++ reader
  AthenaReader *p_athena_reader;
  try
  {
    p_athena_reader = new AthenaReader(p_input_reader);
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
    std::cout << "Error: Could not set up AthenaReader.\n";
    return 1;
  }

  // Set up radiation integrator
  RadiationIntegrator *p_radiation_integrator;
  try
  {
    p_radiation_integrator =
        new RadiationIntegrator(p_input_reader, p_geodesic_integrator, p_athena_reader);
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (const std::bad_optional_access &exception)
  {
    std::cout << "Error: RadiationIntegrator unable to find all needed values in input file.\n";
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not set up RadiationIntegrator.\n";
    return 1;
  }

  // Set up output writer
  OutputWriter *p_output_writer;
  try
  {
    p_output_writer =
        new OutputWriter(p_input_reader, p_geodesic_integrator, p_radiation_integrator);
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
    std::cout << "Error: Could not set up OutputWriter.\n";
    return 1;
  }

  // Go through runs
  for (int n = 0; n < num_runs; n++)
  {
    // Adjust file names
    p_input_reader->AdjustFileNames(n);

    // Read Athena++ file
    try
    {
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
      std::cout << "Error: Could not read Athena++ file.\n";
      return 1;
    }

    // Integrate radiation
    try
    {
      p_radiation_integrator->Integrate();
    }
    catch (const BlacklightException &exception)
    {
      std::cout << exception.what();
      return 1;
    }
    catch (const std::bad_optional_access &exception)
    {
      std::cout << "Error: RadiationIntegrator unable to find all needed values in input file.\n";
      return 1;
    }
    catch (...)
    {
      std::cout << "Error: Could not integrate radiation.\n";
      return 1;
    }

    // Write output
    try
    {
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
  }

  // Free memory
  delete p_output_writer;
  delete p_radiation_integrator;
  delete p_athena_reader;
  delete p_geodesic_integrator;
  delete p_input_reader;

  // End program
  return 0;
}
