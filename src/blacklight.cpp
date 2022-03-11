// Blacklight main file

// C++ headers
#include <iomanip>   // setprecision
#include <iostream>  // cout
#include <optional>  // bad_optional_access, optional
#include <string>    // string

// Library headers
#include <omp.h>  // omp_get_wtime, omp_set_num_threads

// Blacklight headers
#include "blacklight.hpp"
#include "geodesic_integrator/geodesic_integrator.hpp"    // GeodesicIntegrator
#include "input_reader/input_reader.hpp"                  // InputReader
#include "output_writer/output_writer.hpp"                // OutputWriter
#include "radiation_integrator/radiation_integrator.hpp"  // RadiationIntegrator
#include "simulation_reader/simulation_reader.hpp"        // SimulationReader
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
  // Prepare timers
  double time_start = omp_get_wtime();
  double time_geodesic = 0.0;
  double time_read = 0.0;
  double time_sample = 0.0;
  double time_image = 0.0;
  double time_render = 0.0;

  // Parse command-line inputs
  if (argc != 2)
  {
    std::cout << "Error: Must give a single input file.\n";
    return 1;
  }
  const std::string input_file(argv[1]);

  // Prepare pointers to objects
  InputReader *p_input_reader;
  GeodesicIntegrator *p_geodesic_integrator;
  SimulationReader *p_simulation_reader;
  RadiationIntegrator *p_radiation_integrator;
  OutputWriter *p_output_writer;

  // Read input file
  int num_runs;
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
  try
  {
    p_geodesic_integrator = new GeodesicIntegrator(p_input_reader);
    time_geodesic += p_geodesic_integrator->Integrate();
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

  // Set up simulation reader
  try
  {
    p_simulation_reader = new SimulationReader(p_input_reader);
  }
  catch (const BlacklightException &exception)
  {
    std::cout << exception.what();
    return 1;
  }
  catch (const std::bad_optional_access &exception)
  {
    std::cout << "Error: SimulationReader unable to find all needed values in input file.\n";
    return 1;
  }
  catch (...)
  {
    std::cout << "Error: Could not set up SimulationReader.\n";
    return 1;
  }

  // Set up radiation integrator
  try
  {
    p_radiation_integrator =
        new RadiationIntegrator(p_input_reader, p_geodesic_integrator, p_simulation_reader);
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
    // Read simulation file
    try
    {
      time_read += p_simulation_reader->Read(n);
    }
    catch (const BlacklightException &exception)
    {
      std::cout << exception.what();
      return 1;
    }
    catch (...)
    {
      std::cout << "Error: Could not read simulation file.\n";
      return 1;
    }

    // Iterate with adaptive refinement
    bool adaptive_complete = false;
    while (not adaptive_complete)
    {
      // Integrate radiation
      try
      {
        adaptive_complete =
            p_radiation_integrator->Integrate(n, &time_sample, &time_image, &time_render);
      }
      catch (const BlacklightException &exception)
      {
        std::cout << exception.what();
        return 1;
      }
      catch (...)
      {
        std::cout << "Error: Could not integrate radiation.\n";
        return 1;
      }

      // Sample additional geodesics
      if (not adaptive_complete)
        try
        {
          time_geodesic += p_geodesic_integrator->AddGeodesics(p_radiation_integrator);
        }
        catch (const BlacklightException &exception)
        {
          std::cout << exception.what();
          return 1;
        }
        catch (...)
        {
          std::cout << "Error: Could not integrate geodesics.\n";
          return 1;
        }
    }

    // Write output
    try
    {
      p_output_writer->Write(n);
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
  }

  // Free memory
  delete p_output_writer;
  delete p_radiation_integrator;
  delete p_simulation_reader;
  delete p_geodesic_integrator;
  delete p_input_reader;

  // Report timings
  double time_full = omp_get_wtime() - time_start;
  std::cout << std::setprecision(7);
  std::cout << "\nCalculation completed.";
  std::cout << "\nElapsed time:            " << time_full << " s";
  std::cout << "\n  Integrating geodesics: " << time_geodesic << " s";
  std::cout << "\n  Reading simulation:    " << time_read << " s";
  std::cout << "\n  Sampling simulation:   " << time_sample << " s";
  std::cout << "\n  Integrating image:     " << time_image << " s";
  std::cout << "\n  Rendering:             " << time_render << " s";
  std::cout << "\n\n";

  // End program
  return 0;
}
