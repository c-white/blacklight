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
//   Saves certain sample data (sample_inds, sample_fracs (if needed), sample_nan, sample_fallback).
void RadiationIntegrator::SaveSampling()
{
  // Open checkpoint file for writing
  std::ofstream checkpoint_stream(checkpoint_sample_file,
      std::ios_base::out | std::ios_base::binary);
  if (not checkpoint_stream.is_open())
    throw BlacklightException("Could not open sample checkpoint file.");

  // Write sampling data
  WriteBinary(&checkpoint_stream, sample_inds);
  if (simulation_interp)
    WriteBinary(&checkpoint_stream, sample_fracs);
  WriteBinary(&checkpoint_stream, sample_nan);
  WriteBinary(&checkpoint_stream, sample_fallback);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for loading sampling data
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Reads file specified by checkpoint_sample_file.
//   Saves certain sample data (sample_inds, sample_fracs (if needed), sample_nan, sample_fallback),
//       allocating arrays.
void RadiationIntegrator::LoadSampling()
{
  // Open checkpoint file for readiing
  std::ifstream checkpoint_stream(checkpoint_sample_file,
      std::ios_base::in | std::ios_base::binary);
  if (not checkpoint_stream.is_open())
    throw BlacklightException("Could not open sample checkpoint file.");

  // Read sampling data
  ReadBinary(&checkpoint_stream, &sample_inds);
  if (simulation_interp)
    ReadBinary(&checkpoint_stream, &sample_fracs);
  ReadBinary(&checkpoint_stream, &sample_nan);
  ReadBinary(&checkpoint_stream, &sample_fallback);
  return;
}
