// Blacklight geodesic integrator

// C++ headers
#include <cmath>     // abs, acos, cos, sqrt
#include <optional>  // optional
#include <string>    // string

// Library headers
#include <omp.h>  // omp_get_wtime

// Blacklight headers
#include "geodesic_integrator.hpp"
#include "../blacklight.hpp"                                 // enums
#include "../input_reader/input_reader.hpp"                  // InputReader
#include "../radiation_integrator/radiation_integrator.hpp"  // RadiationIntegrator
#include "../utils/exceptions.hpp"                           // BlacklightException

//--------------------------------------------------------------------------------------------------

// Geodesic integrator constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
GeodesicIntegrator::GeodesicIntegrator(const InputReader *p_input_reader)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();

  // Copy checkpoint parameters
  checkpoint_geodesic_save = p_input_reader->checkpoint_geodesic_save.value();
  checkpoint_geodesic_load = p_input_reader->checkpoint_geodesic_load.value();
  if (checkpoint_geodesic_save and checkpoint_geodesic_load)
    throw BlacklightException("Cannot both save and load a geodesic checkpoint.");
  if (checkpoint_geodesic_save or checkpoint_geodesic_load)
    checkpoint_geodesic_file = p_input_reader->checkpoint_geodesic_file.value();

  // Copy camera parameters
  camera_type = p_input_reader->camera_type.value();
  camera_r = p_input_reader->camera_r.value();
  camera_th = p_input_reader->camera_th.value();
  camera_ph = p_input_reader->camera_ph.value();
  camera_urn = p_input_reader->camera_urn.value();
  camera_uthn = p_input_reader->camera_uthn.value();
  camera_uphn = p_input_reader->camera_uphn.value();
  camera_k_r = p_input_reader->camera_k_r.value();
  camera_k_th = p_input_reader->camera_k_th.value();
  camera_k_ph = p_input_reader->camera_k_ph.value();
  camera_rotation = p_input_reader->camera_rotation.value();
  camera_width = p_input_reader->camera_width.value();
  camera_resolution = p_input_reader->camera_resolution.value();
  camera_pole = p_input_reader->camera_pole.value();

  // Copy ray-tracing parameters
  ray_flat = p_input_reader->ray_flat.value();
  ray_terminate = p_input_reader->ray_terminate.value();
  if (ray_terminate != RayTerminate::photon)
    ray_factor = p_input_reader->ray_factor.value();
  ray_step = p_input_reader->ray_step.value();
  ray_max_steps = p_input_reader->ray_max_steps.value();
  ray_max_retries = p_input_reader->ray_max_retries.value();
  ray_tol_abs = p_input_reader->ray_tol_abs.value();
  ray_tol_rel = p_input_reader->ray_tol_rel.value();
  ray_err_factor = p_input_reader->ray_err_factor.value();
  ray_min_factor = p_input_reader->ray_min_factor.value();
  ray_max_factor = p_input_reader->ray_max_factor.value();

  // Copy image parameters
  image_frequency = p_input_reader->image_frequency.value();
  image_normalization = p_input_reader->image_normalization.value();

  // Copy adaptive parameters
  adaptive_max_level = p_input_reader->adaptive_max_level.value();
  if (adaptive_max_level > 0)
  {
    adaptive_block_size = p_input_reader->adaptive_block_size.value();
    if (adaptive_block_size <= 0)
      throw BlacklightException("Must have positive adaptive_block_size.");
    if (camera_resolution % adaptive_block_size != 0)
      throw BlacklightException("Must have adaptive_block_size divide camera_resolution.");
  }

  // Set and calculate geometry data
  if (model_type == ModelType::simulation)
  {
    bh_m = 1.0;
    bh_a = p_input_reader->simulation_a.value();
  }
  else if (model_type == ModelType::formula)
  {
    bh_m = 1.0;
    bh_a = p_input_reader->formula_spin.value();
  }
  double r_horizon = bh_m + std::sqrt(bh_m * bh_m - bh_a * bh_a);
  if (ray_terminate == RayTerminate::photon)
    r_terminate = 2.0 * bh_m * (1.0 + std::cos(2.0 / 3.0 * std::acos(-std::abs(bh_a) / bh_m)));
  else if (ray_terminate == RayTerminate::multiplicative)
    r_terminate = r_horizon * ray_factor;
  else if (ray_terminate == RayTerminate::additive)
    r_terminate = r_horizon + ray_factor;

  // Calculate number of pixels
  camera_num_pix = camera_resolution * camera_resolution;

  // Allocate space for camera data
  camera_loc = new Array<int>[adaptive_max_level+1];
  camera_pos = new Array<double>[adaptive_max_level+1];
  camera_dir = new Array<double>[adaptive_max_level+1];

  // Allocate space for geodesic data
  geodesic_num_steps = new int[adaptive_max_level+1];
  sample_flags = new Array<bool>[adaptive_max_level+1];
  sample_num = new Array<int>[adaptive_max_level+1];
  sample_pos = new Array<double>[adaptive_max_level+1];
  sample_dir = new Array<double>[adaptive_max_level+1];
  sample_len = new Array<double>[adaptive_max_level+1];

  // Prepare bookkeeping for adaptive refinement
  if (adaptive_max_level > 0)
  {
    linear_root_blocks = camera_resolution / adaptive_block_size;
    block_num_pix = adaptive_block_size * adaptive_block_size;
    camera_loc[0].Allocate(linear_root_blocks * linear_root_blocks, 2);
    for (int block_v = 0, block = 0; block_v < linear_root_blocks; block_v++)
      for (int block_u = 0; block_u < linear_root_blocks; block_u++, block++)
      {
        camera_loc[0](block,0) = block_v;
        camera_loc[0](block,1) = block_u;
      }
  }
}

//--------------------------------------------------------------------------------------------------

// Geodesic integrator destructor
GeodesicIntegrator::~GeodesicIntegrator()
{
  for (int level = 0; level <= adaptive_max_level; level++)
  {
    camera_loc[level].Deallocate();
    camera_pos[level].Deallocate();
    camera_dir[level].Deallocate();
    sample_flags[level].Deallocate();
    sample_num[level].Deallocate();
    sample_pos[level].Deallocate();
    sample_dir[level].Deallocate();
    sample_len[level].Deallocate();
  }
  delete[] camera_loc;
  delete[] camera_pos;
  delete[] camera_dir;
  delete[] geodesic_num_steps;
  delete[] sample_flags;
  delete[] sample_num;
  delete[] sample_pos;
  delete[] sample_dir;
  delete[] sample_len;
}

//--------------------------------------------------------------------------------------------------

// Top-level function for integrating geodesics
// Inputs: (none)
// Outputs:
//   returned value: execution time in seconds
double GeodesicIntegrator::Integrate()
{
  // Prepare timer
  double time_start = omp_get_wtime();

  // Reset adaptive level counter
  adaptive_level = 0;

  // Load data from checkpoint
  if (checkpoint_geodesic_load)
    LoadGeodesics();

  // Calculate geodesics
  if (not checkpoint_geodesic_load)
  {
    InitializeCamera();
    InitializeGeodesics();
    IntegrateGeodesics();
    ReverseGeodesics();
  }

  // Save data to checkpoint
  if (checkpoint_geodesic_save)
    SaveGeodesics();

  // Calculate elapsed time
  return omp_get_wtime() - time_start;
}

//--------------------------------------------------------------------------------------------------

// Function for augmenting rays at new refinement level
// Inputs:
//   p_radiation_integrator: pointer to object containing processed image
// Outputs:
//   returned value: execution time in seconds
// Notes:
//   Acquires values from RadiationIntegrator that were not available at construction.
double GeodesicIntegrator::AddGeodesics(const RadiationIntegrator *p_radiation_integrator)
{
  // Prepare timer
  double time_start = omp_get_wtime();

  // Increment adaptive level counter
  adaptive_level++;

  // Acquire refinement data
  block_counts = p_radiation_integrator->block_counts;
  refinement_flags = p_radiation_integrator->refinement_flags;

  // Calculate geodesics
  AugmentCamera();
  InitializeGeodesics();
  IntegrateGeodesics();
  ReverseGeodesics();

  // Calculate elapsed time
  return omp_get_wtime() - time_start;
}
