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
//   Saves camera data (momentum_factor, cam_x, u_con, u_cov, norm_con, norm_con_c, hor_con_c,
//       vert_con_c, camera_pos[0], and camera_dir[0]).
//   Does not save camera_num_pix, which is calculated by constructor.
//   Saves geodesic data (geodesic_num_steps[0], sample_flags[0], sample_num[0], sample_pos[0],
//       sample_dir[0], and sample_len[0]).
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
  WriteBinary(&checkpoint_stream, cam_x, 4);
  WriteBinary(&checkpoint_stream, u_con, 4);
  WriteBinary(&checkpoint_stream, u_cov, 4);
  WriteBinary(&checkpoint_stream, norm_con, 4);
  WriteBinary(&checkpoint_stream, norm_con_c, 4);
  WriteBinary(&checkpoint_stream, hor_con_c, 4);
  WriteBinary(&checkpoint_stream, vert_con_c, 4);
  WriteBinary(&checkpoint_stream, camera_pos[0]);
  WriteBinary(&checkpoint_stream, camera_dir[0]);

  // Write geodesic data
  WriteBinary(&checkpoint_stream, geodesic_num_steps[0]);
  WriteBinary(&checkpoint_stream, sample_flags[0]);
  WriteBinary(&checkpoint_stream, sample_num[0]);
  WriteBinary(&checkpoint_stream, sample_pos[0]);
  WriteBinary(&checkpoint_stream, sample_dir[0]);
  WriteBinary(&checkpoint_stream, sample_len[0]);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for loading geodesic data
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Reads file specified by checkpoint_geodesic_file.
//   Initializes camera data (momentum_factor, cam_x, u_con, u_cov, norm_con, norm_con_c, hor_con_c,
//       vert_con_c, camera_pos[0], and camera_dir[0]), allocating arrays where necessary.
//   Does not initialize camera_num_pix, which is calculated by constructor.
//   Initializes geodesic data (geodesic_num_steps[0], sample_flags[0], sample_num[0],
//       sample_pos[0], sample_dir[0], and sample_len[0]), allocating arrays where necessary.
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
  ReadBinary(&checkpoint_stream, cam_x, 4);
  ReadBinary(&checkpoint_stream, u_con, 4);
  ReadBinary(&checkpoint_stream, u_cov, 4);
  ReadBinary(&checkpoint_stream, norm_con, 4);
  ReadBinary(&checkpoint_stream, norm_con_c, 4);
  ReadBinary(&checkpoint_stream, hor_con_c, 4);
  ReadBinary(&checkpoint_stream, vert_con_c, 4);
  ReadBinary(&checkpoint_stream, &camera_pos[0]);
  ReadBinary(&checkpoint_stream, &camera_dir[0]);

  // Read geodesic data
  ReadBinary(&checkpoint_stream, &geodesic_num_steps[0]);
  ReadBinary(&checkpoint_stream, &sample_flags[0]);
  ReadBinary(&checkpoint_stream, &sample_num[0]);
  ReadBinary(&checkpoint_stream, &sample_pos[0]);
  ReadBinary(&checkpoint_stream, &sample_dir[0]);
  ReadBinary(&checkpoint_stream, &sample_len[0]);
  return;
}
