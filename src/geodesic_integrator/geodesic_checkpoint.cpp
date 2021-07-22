// Blacklight geodesic integrator - checkpoint saving and loading

// C++ headers
#include <fstream>  // ifstream, ofstream
#include <ios>      // ios_base

// Blacklight headers
#include "geodesic_integrator.hpp"
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException
#include "../utils/file_io.hpp"     // ReadBinary, WriteBinary

//--------------------------------------------------------------------------------------------------

// Function for saving geodesic data
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Overwrites file specified by checkpoint_geodesic_file.
//   Saves camera data (momentum_factor, camera_ucon, camera_ucov, camera_up_con_c, camera_num_pix,
//       camera_pos, and camera_dir).
//   Saves geodesic data (geodesic_num_steps, sample_flags, sample_num, sample_pos, sample_dir, and
//       sample_len).
//   Does not save geodesic_pos, geodesic_dir, or geodesic_len, which are deallocated after internal
//       use.
void GeodesicIntegrator::SaveGeodesics()
{
  // Open checkpoint file for writing
  std::ofstream checkpoint_stream(checkpoint_geodesic_file,
      std::ios_base::out | std::ios_base::binary);
  if (not checkpoint_stream.is_open())
    throw BlacklightException("Could not open geodesic checkpoint file.");

  // Write camera data
  WriteBinary(&checkpoint_stream, momentum_factor);
  WriteBinary(&checkpoint_stream, camera_ucon, 4);
  WriteBinary(&checkpoint_stream, camera_ucov, 4);
  WriteBinary(&checkpoint_stream, camera_up_con_c, 4);
  WriteBinary(&checkpoint_stream, camera_pos);
  WriteBinary(&checkpoint_stream, camera_dir);

  // Write geodesic data
  WriteBinary(&checkpoint_stream, geodesic_num_steps);
  WriteBinary(&checkpoint_stream, sample_flags);
  WriteBinary(&checkpoint_stream, sample_num);
  WriteBinary(&checkpoint_stream, sample_pos);
  WriteBinary(&checkpoint_stream, sample_dir);
  WriteBinary(&checkpoint_stream, sample_len);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for loading geodesic data
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Reads file specified by checkpoint_geodesic_file.
//   Initializes camera data (momentum_factor, camera_ucon, camera_ucov, camera_up_con_c,
//       camera_num_pix, camera_pos, and camera_dir), allocating arrays where necessary.
//   Initializes geodesic data (geodesic_num_steps, sample_flags, sample_num, sample_pos,
//       sample_dir, and sample_len), allocating arrays where necessary.
//   Does not initialize geodesic_pos, geodesic_dir, or geodesic_len, which are deallocated after
//       internal use when no checkpoint is being loaded.
void GeodesicIntegrator::LoadGeodesics()
{
  // Open checkpoint file for readiing
  std::ifstream checkpoint_stream(checkpoint_geodesic_file,
      std::ios_base::in | std::ios_base::binary);
  if (not checkpoint_stream.is_open())
    throw BlacklightException("Could not open geodesic checkpoint file.");

  // Read camera data
  ReadBinary(&checkpoint_stream, &momentum_factor);
  ReadBinary(&checkpoint_stream, camera_ucon, 4);
  ReadBinary(&checkpoint_stream, camera_ucov, 4);
  ReadBinary(&checkpoint_stream, camera_up_con_c, 4);
  ReadBinary(&checkpoint_stream, &camera_pos);
  ReadBinary(&checkpoint_stream, &camera_dir);

  // Read geodesic data
  ReadBinary(&checkpoint_stream, &geodesic_num_steps);
  ReadBinary(&checkpoint_stream, &sample_flags);
  ReadBinary(&checkpoint_stream, &sample_num);
  ReadBinary(&checkpoint_stream, &sample_pos);
  ReadBinary(&checkpoint_stream, &sample_dir);
  ReadBinary(&checkpoint_stream, &sample_len);
  return;
}
