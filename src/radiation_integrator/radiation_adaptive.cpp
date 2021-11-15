// Blacklight radiation integrator - adaptive ray tracing

// C++ headers
#include <cmath>  // abs, hypot, isfinite

// Library headers
#include <omp.h>  // pragmas, omp_get_thread_num

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for checking refinement of all blocks at most recently integrated refinement level
// Inputs: (none)
// Outputs:
//   returned value: flag indicating no additional geodesics need to be run for this snapshot
bool RadiationIntegrator::CheckAdaptiveRefinement()
{
  // Handle case where no further refinement can be done
  if (adaptive_level >= adaptive_max_level)
    return true;

  // Ensure there is storage for refinement flags
  if (not refinement_flags[adaptive_level].allocated
      or refinement_flags[adaptive_level].n_tot < block_counts[adaptive_level])
  {
    refinement_flags[adaptive_level].Deallocate();
    refinement_flags[adaptive_level].Allocate(block_counts[adaptive_level]);
  }

  // Calculate size of image in blocks at current level
  int linear_num_blocks = linear_root_blocks;
  for (int n = 1; n <= adaptive_level; n++)
    linear_num_blocks *= 2;

  // Work in parallel
  int num_refined_blocks = 0;
  #pragma omp parallel
  {
    // Determine thread number
    int thread = omp_get_thread_num();

    // Go through blocks at root level
    if (adaptive_level == 0)
    {
      #pragma omp for schedule(static)
      for (int block = 0; block < block_counts[0]; block++)
      {
        // Check forced refinement
        if (adaptive_num_regions > 0)
        {
          double y = ((camera_loc[0](block,0) + 0.5) / linear_num_blocks - 0.5)
              * camera_width;
          double x = ((camera_loc[0](block,1) + 0.5) / linear_num_blocks - 0.5)
              * camera_width;
          refinement_flags[0](block) = false;
          for (int n_r = 0; n_r < adaptive_num_regions; n_r++)
            if (0 < adaptive_region_levels[n_r] and x > adaptive_region_x_min_vals[n_r]
                and x < adaptive_region_x_max_vals[n_r] and y > adaptive_region_y_min_vals[n_r]
                and y < adaptive_region_y_max_vals[n_r])
            {
              refinement_flags[0](block) = true;
              break;
            }
          if (refinement_flags[0](block))
            continue;
        }

        // Check remaining refinement conditions
        int s_full_start = adaptive_frequency_num
            * (model_type == ModelType::simulation and image_polarization ? 4 : 1);
        int s_full_end =
            s_full_start + (model_type == ModelType::simulation and image_polarization ? 4 : 1);
        int j_full_start = block / linear_root_blocks * adaptive_block_size;
        int i_full_start = block % linear_root_blocks * adaptive_block_size;
        for (int s = 0, s_full = s_full_start; s_full < s_full_end; s++, s_full++)
          for (int j = 0, j_full = j_full_start; j < adaptive_block_size; j++, j_full++)
            for (int i = 0, i_full = i_full_start; i < adaptive_block_size; i++, i_full++)
            {
              int m = j_full * camera_resolution + i_full;
              image_blocks[thread](s,j,i) = image[0](s_full,m);
            }
        refinement_flags[0](block) = EvaluateBlock(thread);
      }
    }

    // Go through blocks beyond root level
    else
    {
      #pragma omp for schedule(static)
      for (int block = 0; block < block_counts[adaptive_level]; block++)
      {
        // Check forced refinement
        if (adaptive_num_regions > 0)
        {
          double y = ((camera_loc[adaptive_level](block,0) + 0.5) / linear_num_blocks - 0.5)
              * camera_width;
          double x = ((camera_loc[adaptive_level](block,1) + 0.5) / linear_num_blocks - 0.5)
              * camera_width;
          refinement_flags[adaptive_level](block) = false;
          for (int n_r = 0; n_r < adaptive_num_regions; n_r++)
            if (adaptive_level < adaptive_region_levels[n_r] and x > adaptive_region_x_min_vals[n_r]
                and x < adaptive_region_x_max_vals[n_r] and y > adaptive_region_y_min_vals[n_r]
                and y < adaptive_region_y_max_vals[n_r])
            {
              refinement_flags[adaptive_level](block) = true;
              break;
            }
          if (refinement_flags[adaptive_level](block))
            continue;
        }

        // Check remaining refinement conditions
        int s_full_start = adaptive_frequency_num
            * (model_type == ModelType::simulation and image_polarization ? 4 : 1);
        int s_full_end =
            s_full_start + (model_type == ModelType::simulation and image_polarization ? 4 : 1);
        int m_start = block * block_num_pix;
        for (int s = 0, s_full = s_full_start; s_full < s_full_end; s++, s_full++)
          for (int j = 0, m = m_start; j < adaptive_block_size; j++)
            for (int i = 0; i < adaptive_block_size; i++, m++)
              image_blocks[thread](s,j,i) = image[adaptive_level](s_full,m);
        refinement_flags[adaptive_level](block) = EvaluateBlock(thread);
      }
    }

    // Calculate number of blocks needed for next level
    #pragma omp for schedule(static) reduction(+: num_refined_blocks)
    for (int block = 0; block < block_counts[adaptive_level]; block++)
      if (refinement_flags[adaptive_level](block))
        num_refined_blocks++;
  }

  // Record number of blocks needed for next level
  block_counts[adaptive_level+1] = num_refined_blocks * 4;
  return num_refined_blocks == 0;
}

//--------------------------------------------------------------------------------------------------

// Function for determining if a block needs to be refined
// Inputs:
//   thread: thread number indicating which block should be examined
// Outputs:
//   returned value: flag indicating block needs to be refined
// Notes:
//   There are up to 5 similar evaluations. In each case, a quantity Q is computed on n points or
//       (overlapping) chunks of points. Let k be the number of times Q > C for the user-specified
//       cut C. If k/n > F for the user-specified fraction F, the block is flagged for refinement.
//       If F < 0, the test is not run. If F = 0, any point or chunk with Q > C will trigger
//       refinement. A block will be flagged for refinement if any test that is run triggers
//       refinement.
//   The possible definitions of Q are:
//     (1) |I_nu|,
//     (2) |grad(I_nu)|,
//     (3) |grad(I_nu) / I_nu|,
//     (4) |lapl(I_nu)|,
//     (5) |lapl(I_nu) / I_nu|.
//   Lengths used in the above derivatives are taken to be units of separation between points; that
//       is, options (2)-(4) should decrease as refinement level increases.
bool RadiationIntegrator::EvaluateBlock(int thread)
{
  // Extract intensity
  Array<double> intensity = image_blocks[thread];
  intensity.Slice(3, 0, 0);

  // Test value of intensity
  if (adaptive_val_frac >= 0.0)
  {
    int num_examined = 0;
    int num_exceeded = 0;
    for (int i = 0; i < adaptive_block_size; i++)
      for (int j = 0; j < adaptive_block_size; j++)
      {
        double q = std::abs(intensity(i,j));
        if (std::isfinite(q))
          num_examined++;
        else
          continue;
        if (q > adaptive_val_cut)
          num_exceeded++;
      }
    double frac = static_cast<double>(num_exceeded) / static_cast<double>(num_examined);
    if (frac > adaptive_val_frac)
      return true;
  }

  // Test absolute gradient of intensity
  if (adaptive_abs_grad_frac >= 0.0)
  {
    int num_examined = 0;
    int num_exceeded = 0;
    for (int i = 0; i < adaptive_block_size; i++)
      for (int j = 0; j < adaptive_block_size; j++)
      {
        double q_x;
        if (j == 0)
          q_x = intensity(i,j+1) - intensity(i,j);
        else if (j == adaptive_block_size - 1)
          q_x = intensity(i,j) - intensity(i,j-1);
        else
          q_x = 0.5 * (intensity(i,j+1) - intensity(i,j-1));
        double q_y;
        if (i == 0)
          q_y = intensity(i+1,j) - intensity(i,j);
        else if (i == adaptive_block_size - 1)
          q_y = intensity(i,j) - intensity(i-1,j);
        else
          q_y = 0.5 * (intensity(i+1,j) - intensity(i-1,j));
        double q = std::hypot(q_x, q_y);
        if (std::isfinite(q))
          num_examined++;
        else
          continue;
        if (q > adaptive_abs_grad_cut)
          num_exceeded++;
      }
    double frac = static_cast<double>(num_exceeded) / static_cast<double>(num_examined);
    if (frac > adaptive_abs_grad_frac)
      return true;
  }

  // Test relative gradient of intensity
  if (adaptive_rel_grad_frac >= 0.0)
  {
    int num_examined = 0;
    int num_exceeded = 0;
    for (int i = 0; i < adaptive_block_size; i++)
      for (int j = 0; j < adaptive_block_size; j++)
      {
        double q_x;
        if (j == 0)
          q_x = 2.0 * (intensity(i,j+1) - intensity(i,j)) / (intensity(i,j) + intensity(i,j+1));
        else if (j == adaptive_block_size - 1)
          q_x = 2.0 * (intensity(i,j) - intensity(i,j-1)) / (intensity(i,j-1) + intensity(i,j));
        else
          q_x = 2.0 * (intensity(i,j+1) - intensity(i,j-1))
              / (intensity(i,j-1) + 2.0 * intensity(i,j) + intensity(i,j+1));
        double q_y;
        if (i == 0)
          q_y = 2.0 * (intensity(i+1,j) - intensity(i,j)) / (intensity(i,j) + intensity(i+1,j));
        else if (i == adaptive_block_size - 1)
          q_y = 2.0 * (intensity(i,j) - intensity(i-1,j)) / (intensity(i-1,j) + intensity(i,j));
        else
          q_y = 2.0 * (intensity(i+1,j) - intensity(i-1,j))
              / (intensity(i-1,j) + 2.0 * intensity(i,j) + intensity(i+1,j));
        double q = std::hypot(q_x, q_y);
        if (std::isfinite(q))
          num_examined++;
        else
          continue;
        if (q > adaptive_rel_grad_cut)
          num_exceeded++;
      }
    double frac = static_cast<double>(num_exceeded) / static_cast<double>(num_examined);
    if (frac > adaptive_rel_grad_frac)
      return true;
  }

  // Test absolute Laplacian of intensity
  if (adaptive_abs_lapl_frac >= 0.0)
  {
    int num_examined = 0;
    int num_exceeded = 0;
    for (int i = 1; i < adaptive_block_size - 1; i++)
      for (int j = 1; j < adaptive_block_size - 1; j++)
      {
        double q_x = intensity(i,j-1) - 2.0 * intensity(i,j) + intensity(i,j+1);
        double q_y = intensity(i-1,j) - 2.0 * intensity(i,j) + intensity(i+1,j);
        double q = std::abs(q_x + q_y);
        if (std::isfinite(q))
          num_examined++;
        else
          continue;
        if (q > adaptive_abs_lapl_cut)
          num_exceeded++;
      }
    double frac = static_cast<double>(num_exceeded) / static_cast<double>(num_examined);
    if (frac > adaptive_abs_lapl_frac)
      return true;
  }

  // Test relative Laplacian of intensity
  if (adaptive_rel_lapl_frac >= 0.0)
  {
    int num_examined = 0;
    int num_exceeded = 0;
    for (int i = 1; i < adaptive_block_size - 1; i++)
      for (int j = 1; j < adaptive_block_size - 1; j++)
      {
        double q_x = 4.0 * (intensity(i,j-1) - 2.0 * intensity(i,j) + intensity(i,j+1))
            / (intensity(i,j-1) + 2.0 * intensity(i,j) + intensity(i,j+1));
        double q_y = 4.0 * (intensity(i-1,j) - 2.0 * intensity(i,j) + intensity(i+1,j))
            / (intensity(i-1,j) + 2.0 * intensity(i,j) + intensity(i+1,j));
        double q = std::abs(q_x + q_y);
        if (std::isfinite(q))
          num_examined++;
        else
          continue;
        if (q > adaptive_rel_lapl_cut)
          num_exceeded++;
      }
    double frac = static_cast<double>(num_exceeded) / static_cast<double>(num_examined);
    if (frac > adaptive_rel_lapl_frac)
      return true;
  }

  // Conclude no refinement is needed
  return false;
}
