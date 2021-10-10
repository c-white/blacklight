// Blacklight radiation integrator - checkpoint saving and loading

// C++ headers
#include <fstream>  // ifstream, ofstream
#include <ios>      // ios_base

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../utils/array.hpp"        // Array
#include "../utils/exceptions.hpp"   // BlacklightException
#include "../utils/file_io.hpp"      // ReadBinary, WriteBinary

//--------------------------------------------------------------------------------------------------

// Function for saving sampling data
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Overwrites file specified by checkpoint_sample_file.
//   Saves certain sample data (sample_inds[0], sample_fracs[0] (if needed), sample_nan[0], and
//       sample_fallback[0]).
void RadiationIntegrator::SaveSampling()
{
  // Open checkpoint file for writing
  std::ofstream checkpoint_stream(checkpoint_sample_file,
      std::ios_base::out | std::ios_base::binary);
  if (not checkpoint_stream.is_open())
    throw BlacklightException("Could not open sample checkpoint file.");

  // Write sampling data
  WriteBinary(&checkpoint_stream, sample_inds[0]);
  if (simulation_interp)
    WriteBinary(&checkpoint_stream, sample_fracs[0]);
  WriteBinary(&checkpoint_stream, sample_nan[0]);
  WriteBinary(&checkpoint_stream, sample_fallback[0]);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for loading sampling data
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Reads file specified by checkpoint_sample_file.
//   Saves certain sample data (sample_inds[0], sample_fracs[0] (if needed), sample_nan[0], and
//       sample_fallback[0]), allocating arrays.
void RadiationIntegrator::LoadSampling()
{
  // Open checkpoint file for readiing
  std::ifstream checkpoint_stream(checkpoint_sample_file,
      std::ios_base::in | std::ios_base::binary);
  if (not checkpoint_stream.is_open())
    throw BlacklightException("Could not open sample checkpoint file.");

  // Read sampling data
  ReadBinary(&checkpoint_stream, &sample_inds[0]);
  if (simulation_interp)
    ReadBinary(&checkpoint_stream, &sample_fracs[0]);
  ReadBinary(&checkpoint_stream, &sample_nan[0]);
  ReadBinary(&checkpoint_stream, &sample_fallback[0]);
  return;
}
