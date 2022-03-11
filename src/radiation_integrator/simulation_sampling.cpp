// Blacklight radiation integrator - simulation sampling

// C++ headers
#include <algorithm>  // max, min
#include <limits>     // numeric_limits
#include <sstream>    // ostringstream

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"                           // enums
#include "../simulation_reader/simulation_reader.hpp"  // SimulationReader
#include "../utils/array.hpp"                          // Array
#include "../utils/exceptions.hpp"                     // BlacklightException, BlacklightWarning

//--------------------------------------------------------------------------------------------------

// Function for making shallow copies of arrays containing simulation data.
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Acquires values from SimulationReader that were not available at construction.
void RadiationIntegrator::ObtainGridData()
{
  // Copy grid metadata
  if (simulation_coord == Coordinates::sph_ks and simulation_interp and simulation_block_interp)
    n_3_root = p_simulation_reader->n_3_root;

  // Copy grid layout
  if (simulation_interp and simulation_block_interp)
  {
    levels = p_simulation_reader->levels;
    locations = p_simulation_reader->locations;
  }

  // Copy coordinates
  x1f = p_simulation_reader->x1f;
  x2f = p_simulation_reader->x2f;
  x3f = p_simulation_reader->x3f;
  x1v = p_simulation_reader->x1v;
  x2v = p_simulation_reader->x2v;
  x3v = p_simulation_reader->x3v;

  // Copy time
  time = p_simulation_reader->time;

  // Copy cell values
  grid_prim = p_simulation_reader->prim;
  grid_bb = p_simulation_reader->bb;

  // Copy indices
  ind_rho = p_simulation_reader->ind_rho;
  ind_pgas = p_simulation_reader->ind_pgas;
  ind_kappa = p_simulation_reader->ind_kappa;
  ind_uu1 = p_simulation_reader->ind_uu1;
  ind_uu2 = p_simulation_reader->ind_uu2;
  ind_uu3 = p_simulation_reader->ind_uu3;
  ind_bb1 = p_simulation_reader->ind_bb1;
  ind_bb2 = p_simulation_reader->ind_bb2;
  ind_bb3 = p_simulation_reader->ind_bb3;

  // Calculate maximum refinement level and number of blocks in x^3-direction at each level
  if (simulation_coord == Coordinates::sph_ks and simulation_interp and simulation_block_interp)
  {
    int n_b = x1f.n2;
    int n_k = x3v.n1;
    max_level = 0;
    for (int b = 0; b < n_b; b++)
      max_level = std::max(max_level, levels(b));
    n_3_level.Allocate(max_level + 1);
    n_3_level(0) = n_3_root / n_k;
    for (int level = 1; level <= max_level; level++)
      n_3_level(level) = n_3_level(level-1) * 2;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for determining how to sample cell data onto rays.
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be processed
// Outputs: (none)
// Notes:
//   Assumes geodesic_num_steps[adaptive_level], sample_flags[adaptive_level],
//       sample_num[adaptive_level], and sample_pos[adaptive_level] have been set.
//   Allocates and initializes sample_inds[adaptive_level], sample_nan[adaptive_level],
//       sample_cut[adaptive_level], and sample_fallback[adaptive_level].
//   Allocates and initializes sample_fracs[adaptive_level] if simulation_interp == true or if
//       slow_light_on == true and slow_interp == true.
//   If simulation_interp == false, locates cell containing geodesic sample point.
//   If simulation_interp == true and simulation_block_interp == false, prepares trilinear
//       interpolation to geodesic sample point from cell centers, using only data within the same
//       block of cells (i.e. sometimes using extrapolation near block edges).
//   If simulation_interp == true and simulation_block_interp == true, prepares trilinear
//       interpolation after obtaining anchor points possibly from neighboring blocks, even at
//       different refinement levels, or across the periodic boundary in spherical coordinates.
//   If slow_light_on == true, chooses from multiple available time slices for each point.
//   If slow_light_on == true and slow_interp == true, prepares interpolation between adjacent (or
//       sometimes identical) time slices.
void RadiationIntegrator::CalculateSimulationSampling(int snapshot)
{
  // Calculate time of snapshot
  double snapshot_time = 0.0;
  if (slow_light_on)
    snapshot_time = slow_t_start + slow_dt * snapshot;

  // Allocate arrays
  int num_pix = camera_num_pix;
  int num_interp_inds = 4;
  if (slow_light_on)
    num_interp_inds++;
  int num_interp_fracs = 0;
  if (simulation_interp)
    num_interp_fracs += 3;
  if (slow_light_on and slow_interp)
    num_interp_fracs++;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
  {
    if (simulation_interp and simulation_block_interp)
      sample_inds[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level], 8,
          num_interp_inds);
    else
      sample_inds[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level],
          num_interp_inds);
    if (num_interp_fracs > 0)
      sample_fracs[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level],
          num_interp_fracs);
    sample_nan[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_cut[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_fallback[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
  }
  sample_nan[adaptive_level].Zero();
  sample_cut[adaptive_level].Zero();
  sample_fallback[adaptive_level].Zero();

  // Prepare bookkeeping for warnings and errors
  int num_extrap_camera_small = 0;
  int num_extrap_camera_large = 0;
  int num_extrap_source_small = 0;
  int num_extrap_source_large = 0;
  double val_extrap_camera_small = 0.0;
  double val_extrap_camera_large = 0.0;
  double val_extrap_source_small = 0.0;
  double val_extrap_source_large = 0.0;

  // Work in parallel
  #pragma omp parallel
  {
    // Prepare bookkeeping
    int n_b = x1f.n2;
    int n_i = x1v.n1;
    int n_j = x2v.n1;
    int n_k = x3v.n1;
    int b = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    double x1_min_block = x1f(b,0);
    double x1_max_block = x1f(b,n_i);
    double x2_min_block = x2f(b,0);
    double x2_max_block = x2f(b,n_j);
    double x3_min_block = x3f(b,0);
    double x3_max_block = x3f(b,n_k);

    // Resample cell data onto geodesics
    #pragma omp for schedule(static) reduction(+: num_extrap_camera_small, \
        num_extrap_camera_large, num_extrap_source_small, num_extrap_source_large) reduction(max: \
        val_extrap_camera_small, val_extrap_camera_large, val_extrap_source_small, \
        val_extrap_source_large)
    for (int m = 0; m < num_pix; m++)
    {
      // Extract number of steps along this geodesic
      int num_steps = sample_num[adaptive_level](m);

      // Set NaN fallback values if geodesic poorly terminated
      if (fallback_nan and sample_flags[adaptive_level](m))
      {
        for (int n = 0; n < num_steps; n++)
          sample_nan[adaptive_level](m,n) = true;
        continue;
      }

      // Prepare bookkeeping for extrapolation
      bool extrap_camera_small = false;
      bool extrap_camera_large = false;
      bool extrap_source_small = false;
      bool extrap_source_large = false;
      double val_extrap_camera_small_local = 0.0;
      double val_extrap_camera_large_local = 0.0;
      double val_extrap_source_small_local = 0.0;
      double val_extrap_source_large_local = 0.0;

      // Go along geodesic
      for (int n = 0; n < num_steps; n++)
      {
        // Extract coordinates
        double x0 = sample_pos[adaptive_level](m,n,0) + snapshot_time;
        double x1 = sample_pos[adaptive_level](m,n,1);
        double x2 = sample_pos[adaptive_level](m,n,2);
        double x3 = sample_pos[adaptive_level](m,n,3);

        // Cut outside camera radius
        double r = RadialGeodesicCoordinate(x1, x2, x3);
        if (r > camera_r)
        {
          sample_cut[adaptive_level](m,n) = true;
          continue;
        }

        // Cut camera plane
        if (cut_omit_near or cut_omit_far)
        {
          double dot_product = x1 * camera_x[1] + x2 * camera_x[2] + x3 * camera_x[3];
          if ((cut_omit_near and dot_product > 0.0) or (cut_omit_far and dot_product < 0.0))
          {
            sample_cut[adaptive_level](m,n) = true;
            continue;
          }
        }

        // Cut spheres
        if ((cut_omit_in >= 0.0 and r < cut_omit_in) or (cut_omit_out >= 0.0 and r > cut_omit_out))
        {
          sample_cut[adaptive_level](m,n) = true;
          continue;
        }

        // Cut arbitrary plane
        if (cut_plane)
        {
          double dot_product = (x1 - cut_plane_origin_x) * cut_plane_normal_x
              + (x2 - cut_plane_origin_y) * cut_plane_normal_y
              + (x3 - cut_plane_origin_z) * cut_plane_normal_z;
          if (dot_product < 0.0)
          {
            sample_cut[adaptive_level](m,n) = true;
            continue;
          }
        }

        // Convert coordinates
        if (simulation_coord == Coordinates::sph_ks)
          CKSToSKS(&x1, &x2, &x3);

        // Calculate time interpolation
        int t_ind = 0;
        double t_frac = 0.0;
        if (slow_light_on)
        {
          if (x0 >= static_cast<double>(time[0]))
          {
            if (x0 > static_cast<double>(time[0] + extrapolation_tolerance))
            {
              extrap_camera_large = true;
              val_extrap_camera_large_local =
                  std::max(val_extrap_camera_large_local, x0 - static_cast<double>(time[0]));
            }
            else if (x0 > static_cast<double>(time[0]))
            {
              extrap_camera_small = true;
              val_extrap_camera_small_local =
                  std::max(val_extrap_camera_small_local, x0 - static_cast<double>(time[0]));
            }
          }
          else if (x0 <= static_cast<double>(time[slow_chunk_size-1]))
          {
            if (x0 < static_cast<double>(time[slow_chunk_size-1] - extrapolation_tolerance))
            {
              extrap_source_large = true;
              val_extrap_source_large_local = std::max(val_extrap_source_large_local,
                  static_cast<double>(time[slow_chunk_size-1]) - x0);
            }
            else if (x0 < static_cast<double>(time[slow_chunk_size-1]))
            {
              extrap_source_small = true;
              val_extrap_source_small_local = std::max(val_extrap_source_small_local,
                  static_cast<double>(time[slow_chunk_size-1]) - x0);
            }
            if (slow_interp)
            {
              t_ind = slow_chunk_size - 2;
              t_frac = 1.0;
            }
            else
              t_ind = slow_chunk_size - 1;
          }
          else
          {
            while (static_cast<double>(time[t_ind++]) > x0);
            t_ind--;
            if (slow_interp)
            {
              t_ind--;
              t_frac = (x0 - static_cast<double>(time[t_ind]))
                  / static_cast<double>(time[t_ind+1] - time[t_ind]);
            }
            else if
                (static_cast<double>(time[t_ind-1]) - x0 <= x0 - static_cast<double>(time[t_ind]))
              t_ind--;
          }
        }

        // Determine block
        if (x1 < x1_min_block or x1 > x1_max_block or x2 < x2_min_block or x2 > x2_max_block
            or x3 < x3_min_block or x3 > x3_max_block)
        {
          // Check if block contains position
          int b_new;
          double x1_min_temp = x1_min_block;
          double x1_max_temp = x1_max_block;
          double x2_min_temp = x2_min_block;
          double x2_max_temp = x2_max_block;
          double x3_min_temp = x3_min_block;
          double x3_max_temp = x3_max_block;
          for (b_new = 0; b_new < n_b; b_new++)
          {
            x1_min_temp = x1f(b_new,0);
            x1_max_temp = x1f(b_new,n_i);
            x2_min_temp = x2f(b_new,0);
            x2_max_temp = x2f(b_new,n_j);
            x3_min_temp = x3f(b_new,0);
            x3_max_temp = x3f(b_new,n_k);
            if (x1 >= x1_min_temp and x1 <= x1_max_temp and x2 >= x2_min_temp and x2 <= x2_max_temp
                and x3 >= x3_min_temp and x3 <= x3_max_temp)
              break;
          }

          // Set fallback values if off grid
          if (b_new == n_b)
          {
            if (fallback_nan)
              sample_nan[adaptive_level](m,n) = true;
            else
              sample_fallback[adaptive_level](m,n) = true;
            continue;
          }

          // Set newly found block as one to search
          b = b_new;
          x1_min_block = x1_min_temp;
          x1_max_block = x1_max_temp;
          x2_min_block = x2_min_temp;
          x2_max_block = x2_max_temp;
          x3_min_block = x3_min_temp;
          x3_max_block = x3_max_temp;
        }

        // Determine cell
        for (i = 0; i < n_i; i++)
          if (static_cast<double>(x1f(b,i+1)) >= x1)
            break;
        for (j = 0; j < n_j; j++)
          if (static_cast<double>(x2f(b,j+1)) >= x2)
            break;
        for (k = 0; k < n_k; k++)
          if (static_cast<double>(x3f(b,k+1)) >= x3)
            break;

        // Prepare to sample values without interpolation
        if (not simulation_interp)
        {
          sample_inds[adaptive_level](m,n,0) = b;
          sample_inds[adaptive_level](m,n,1) = k;
          sample_inds[adaptive_level](m,n,2) = j;
          sample_inds[adaptive_level](m,n,3) = i;
          if (slow_light_on)
            sample_inds[adaptive_level](m,n,4) = t_ind;
          if (slow_light_on and slow_interp)
            sample_fracs[adaptive_level](m,n,0) = t_frac;
        }

        // Prepare to sample values with intrablock interpolation
        else if (not simulation_block_interp)
        {
          int i_m = i == 0 or (i != n_i - 1 and x1 >= static_cast<double>(x1v(b,i))) ? i : i - 1;
          int j_m = j == 0 or (j != n_j - 1 and x2 >= static_cast<double>(x2v(b,j))) ? j : j - 1;
          int k_m = k == 0 or (k != n_k - 1 and x3 >= static_cast<double>(x3v(b,k))) ? k : k - 1;
          double f_i = (x1 - static_cast<double>(x1v(b,i_m)))
              / (static_cast<double>(x1v(b,i_m+1)) - static_cast<double>(x1v(b,i_m)));
          double f_j = (x2 - static_cast<double>(x2v(b,j_m)))
              / (static_cast<double>(x2v(b,j_m+1)) - static_cast<double>(x2v(b,j_m)));
          double f_k = (x3 - static_cast<double>(x3v(b,k_m)))
              / (static_cast<double>(x3v(b,k_m+1)) - static_cast<double>(x3v(b,k_m)));
          sample_inds[adaptive_level](m,n,0) = b;
          sample_inds[adaptive_level](m,n,1) = k_m;
          sample_inds[adaptive_level](m,n,2) = j_m;
          sample_inds[adaptive_level](m,n,3) = i_m;
          if (slow_light_on)
            sample_inds[adaptive_level](m,n,4) = t_ind;
          sample_fracs[adaptive_level](m,n,0) = f_k;
          sample_fracs[adaptive_level](m,n,1) = f_j;
          sample_fracs[adaptive_level](m,n,2) = f_i;
          if (slow_light_on and slow_interp)
            sample_fracs[adaptive_level](m,n,3) = t_frac;
        }

        // Prepare to sample values with interblock interpolation
        else
        {
          // Determine indices to use for interpolation
          int i_m = x1 >= static_cast<double>(x1v(b,i)) ? i : i - 1;
          int j_m = x2 >= static_cast<double>(x2v(b,j)) ? j : j - 1;
          int k_m = x3 >= static_cast<double>(x3v(b,k)) ? k : k - 1;
          int i_p = i_m + 1;
          int j_p = j_m + 1;
          int k_p = k_m + 1;

          // Calculate fractions to use in interpolation
          double x1_m = i_m == -1 ? 2.0 * static_cast<double>(x1f(b,i))
              - static_cast<double>(x1v(b,i)) : static_cast<double>(x1v(b,i_m));
          double x2_m = j_m == -1 ? 2.0 * static_cast<double>(x2f(b,j))
              - static_cast<double>(x2v(b,j)) : static_cast<double>(x2v(b,j_m));
          double x3_m = k_m == -1 ? 2.0 * static_cast<double>(x3f(b,k))
              - static_cast<double>(x3v(b,k)) : static_cast<double>(x3v(b,k_m));
          double x1_p = i_p == n_i ? 2.0 * static_cast<double>(x1v(b,i+1))
              - static_cast<double>(x1v(b,i)) : static_cast<double>(x1v(b,i_p));
          double x2_p = j_p == n_j ? 2.0 * static_cast<double>(x2v(b,j+1))
              - static_cast<double>(x2v(b,j)) : static_cast<double>(x2v(b,j_p));
          double x3_p = k_p == n_k ? 2.0 * static_cast<double>(x3v(b,k+1))
              - static_cast<double>(x3v(b,k)) : static_cast<double>(x3v(b,k_p));
          double f_i = (x1 - x1_m) / (x1_p - x1_m);
          double f_j = (x2 - x2_m) / (x2_p - x2_m);
          double f_k = (x3 - x3_m) / (x3_p - x3_m);

          // Find interpolation anchors
          int inds[8][4];
          FindNearbyInds(b, k_m, j_m, i_m, k, j, i, x3, x2, x1, inds[0]);
          FindNearbyInds(b, k_m, j_m, i_p, k, j, i, x3, x2, x1, inds[1]);
          FindNearbyInds(b, k_m, j_p, i_m, k, j, i, x3, x2, x1, inds[2]);
          FindNearbyInds(b, k_m, j_p, i_p, k, j, i, x3, x2, x1, inds[3]);
          FindNearbyInds(b, k_p, j_m, i_m, k, j, i, x3, x2, x1, inds[4]);
          FindNearbyInds(b, k_p, j_m, i_p, k, j, i, x3, x2, x1, inds[5]);
          FindNearbyInds(b, k_p, j_p, i_m, k, j, i, x3, x2, x1, inds[6]);
          FindNearbyInds(b, k_p, j_p, i_p, k, j, i, x3, x2, x1, inds[7]);

          // Store results
          for (int p = 0; p < 8; p++)
          {
            for (int q = 0; q < 4; q++)
              sample_inds[adaptive_level](m,n,p,q) = inds[p][q];
            if (slow_light_on)
              sample_inds[adaptive_level](m,n,p,4) = t_ind;
          }
          sample_fracs[adaptive_level](m,n,0) = f_k;
          sample_fracs[adaptive_level](m,n,1) = f_j;
          sample_fracs[adaptive_level](m,n,2) = f_i;
          if (slow_light_on and slow_interp)
            sample_fracs[adaptive_level](m,n,3) = t_frac;
        }
      }

      // Add to accounting of warnings and errors
      if (extrap_camera_small)
      {
        num_extrap_camera_small++;
        val_extrap_camera_small = std::max(val_extrap_camera_small, val_extrap_camera_small_local);
      }
      if (extrap_camera_large)
      {
        num_extrap_camera_large++;
        val_extrap_camera_large = std::max(val_extrap_camera_large, val_extrap_camera_large_local);
      }
      if (extrap_source_small)
      {
        num_extrap_source_small++;
        val_extrap_source_small = std::max(val_extrap_source_small, val_extrap_source_small_local);
      }
      if (extrap_source_large)
      {
        num_extrap_source_large++;
        val_extrap_source_large = std::max(val_extrap_source_large, val_extrap_source_large_local);
      }
    }
  }

  // Throw error if large extrapolation needed
  if (num_extrap_camera_large > 0)
  {
    std::ostringstream message;
    message << "Snapshot " << snapshot << " at time " << snapshot_time;
    message << " requires significant extrapolation forward in time (" << num_extrap_camera_large;
    message << "/" << num_pix << " pixels, by up to " << val_extrap_camera_large;
    message << " gravitational times).";
    throw BlacklightException(message.str().c_str());
  }
  if (num_extrap_source_large > 0)
  {
    std::ostringstream message;
    message << "Snapshot " << snapshot << " at time " << snapshot_time;
    message << " requires significant extrapolation backward in time (" << num_extrap_source_large;
    message << "/" << num_pix << " pixels, by up to " << val_extrap_source_large;
    message << " gravitational times).";
    throw BlacklightException(message.str().c_str());
  }

  // Warn if small extrapolation needed
  if (num_extrap_camera_small > 0)
  {
    std::ostringstream message;
    message << "Snapshot " << snapshot << " at time " << snapshot_time;
    message << " requires moderate extrapolation forward in time (" << num_extrap_camera_small;
    message << "/" << num_pix << " pixels, by up to " << val_extrap_camera_small;
    message << " gravitational times).";
    BlacklightWarning(message.str().c_str());
  }
  if (num_extrap_source_small > 0)
  {
    std::ostringstream message;
    message << "Snapshot " << snapshot << " at time " << snapshot_time;
    message << " requires moderate extrapolation backward in time (" << num_extrap_source_small;
    message << "/" << num_pix << " pixels, by up to " << val_extrap_source_small;
    message << " gravitational times).";
    BlacklightWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for resampling simulation cell data onto rays.
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes geodesic_num_steps[adaptive_level], sample_num[adaptive_level],
//       sample_inds[adaptive_level], sample_nan[adaptive_level], sample_cut[adaptive_level],
//       and sample_fallback[adaptive_level] have been set.
//   Assumes sample_fracs[adaptive_level] has been set if simulation_interp == true.
//   Allocates and initializes sample_rho[adaptive_level], sample_pgas[adaptive_level],
//       sample_kappa[adaptive_level] (if needed), sample_uu1[adaptive_level],
//       sample_uu2[adaptive_level], sample_uu3[adaptive_level], sample_bb1[adaptive_level],
//       sample_bb2[adaptive_level], and sample_bb3[adaptive_level].
//   Deallocates sample_inds[adaptive_level], sample_fracs[adaptive_level],
//       sample_nan[adaptive_level], and sample_fallback[adaptive_level] if adaptive_level > 0.
void RadiationIntegrator::SampleSimulation()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
  {
    sample_rho[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_pgas[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    if (plasma_model == PlasmaModel::code_kappa)
      sample_kappa[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_uu1[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_uu2[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_uu3[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_bb1[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_bb2[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
    sample_bb3[adaptive_level].Allocate(num_pix, geodesic_num_steps[adaptive_level]);
  }
  sample_rho[adaptive_level].Zero();
  sample_pgas[adaptive_level].Zero();
  sample_kappa[adaptive_level].Zero();
  sample_uu1[adaptive_level].Zero();
  sample_uu2[adaptive_level].Zero();
  sample_uu3[adaptive_level].Zero();
  sample_bb1[adaptive_level].Zero();
  sample_bb2[adaptive_level].Zero();
  sample_bb3[adaptive_level].Zero();

  // Resample cell data onto geodesics in parallel
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < num_pix; m++)
  {
    // Extract number of steps along this geodesic
    int num_steps = sample_num[adaptive_level](m);

    // Go along geodesic
    for (int n = 0; n < num_steps; n++)
    {
      // Set NaN values
      if (sample_nan[adaptive_level](m,n))
      {
        sample_rho[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_pgas[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        if (plasma_model == PlasmaModel::code_kappa)
          sample_kappa[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_uu1[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_uu2[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_uu3[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_bb1[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_bb2[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
        sample_bb3[adaptive_level](m,n) = std::numeric_limits<float>::quiet_NaN();
      }

      // Skip cut regions
      else if (sample_cut[adaptive_level](m,n))
        continue;

      // Set fallback values
      else if (sample_fallback[adaptive_level](m,n))
      {
        sample_rho[adaptive_level](m,n) = fallback_rho;
        sample_pgas[adaptive_level](m,n) = fallback_pgas;
        if (plasma_model == PlasmaModel::code_kappa)
          sample_kappa[adaptive_level](m,n) = fallback_kappa;
        sample_uu1[adaptive_level](m,n) = fallback_uu1;
        sample_uu2[adaptive_level](m,n) = fallback_uu2;
        sample_uu3[adaptive_level](m,n) = fallback_uu3;
        sample_bb1[adaptive_level](m,n) = fallback_bb1;
        sample_bb2[adaptive_level](m,n) = fallback_bb2;
        sample_bb3[adaptive_level](m,n) = fallback_bb3;
      }

      // Set nearest values
      else if (not simulation_interp)
      {
        // Extract indices
        int b = sample_inds[adaptive_level](m,n,0);
        int k = sample_inds[adaptive_level](m,n,1);
        int j = sample_inds[adaptive_level](m,n,2);
        int i = sample_inds[adaptive_level](m,n,3);
        int t = 0;
        if (slow_light_on)
          t = sample_inds[adaptive_level](m,n,4);

        // Calculate values without temporal interpolation
        if (not (slow_light_on and slow_interp))
        {
          sample_rho[adaptive_level](m,n) = grid_prim[t](ind_rho,b,k,j,i);
          sample_pgas[adaptive_level](m,n) = grid_prim[t](ind_pgas,b,k,j,i);
          if (plasma_model == PlasmaModel::code_kappa)
            sample_kappa[adaptive_level](m,n) = grid_prim[t](ind_kappa,b,k,j,i);
          sample_uu1[adaptive_level](m,n) = grid_prim[t](ind_uu1,b,k,j,i);
          sample_uu2[adaptive_level](m,n) = grid_prim[t](ind_uu2,b,k,j,i);
          sample_uu3[adaptive_level](m,n) = grid_prim[t](ind_uu3,b,k,j,i);
          sample_bb1[adaptive_level](m,n) = grid_bb[t](ind_bb1,b,k,j,i);
          sample_bb2[adaptive_level](m,n) = grid_bb[t](ind_bb2,b,k,j,i);
          sample_bb3[adaptive_level](m,n) = grid_bb[t](ind_bb3,b,k,j,i);
        }

        // Calculate values with temporal interpolation
        else
        {
          // Perform spatial interpolation on first slice
          double rho_1 = static_cast<double>(grid_prim[t](ind_rho,b,k,j,i));
          double pgas_1 = static_cast<double>(grid_prim[t](ind_pgas,b,k,j,i));
          double kappa_1 = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa_1 = static_cast<double>(grid_prim[t](ind_kappa,b,k,j,i));
          double uu1_1 = static_cast<double>(grid_prim[t](ind_uu1,b,k,j,i));
          double uu2_1 = static_cast<double>(grid_prim[t](ind_uu2,b,k,j,i));
          double uu3_1 = static_cast<double>(grid_prim[t](ind_uu3,b,k,j,i));
          double bb1_1 = static_cast<double>(grid_bb[t](ind_bb1,b,k,j,i));
          double bb2_1 = static_cast<double>(grid_bb[t](ind_bb2,b,k,j,i));
          double bb3_1 = static_cast<double>(grid_bb[t](ind_bb3,b,k,j,i));

          // Perform spatial interpolation on second slice
          double rho_2 = static_cast<double>(grid_prim[t+1](ind_rho,b,k,j,i));
          double pgas_2 = static_cast<double>(grid_prim[t+1](ind_pgas,b,k,j,i));
          double kappa_2 = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa_2 = static_cast<double>(grid_prim[t+1](ind_kappa,b,k,j,i));
          double uu1_2 = static_cast<double>(grid_prim[t+1](ind_uu1,b,k,j,i));
          double uu2_2 = static_cast<double>(grid_prim[t+1](ind_uu2,b,k,j,i));
          double uu3_2 = static_cast<double>(grid_prim[t+1](ind_uu3,b,k,j,i));
          double bb1_2 = static_cast<double>(grid_bb[t+1](ind_bb1,b,k,j,i));
          double bb2_2 = static_cast<double>(grid_bb[t+1](ind_bb2,b,k,j,i));
          double bb3_2 = static_cast<double>(grid_bb[t+1](ind_bb3,b,k,j,i));

          // Assign interpolated values
          double t_frac = sample_fracs[adaptive_level](m,n,0);
          sample_rho[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * rho_1 + t_frac * rho_2);
          sample_pgas[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * pgas_1 + t_frac * pgas_2);
          if (plasma_model == PlasmaModel::code_kappa)
            sample_kappa[adaptive_level](m,n) =
                static_cast<float>((1.0 - t_frac) * kappa_1 + t_frac * kappa_2);
          sample_uu1[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu1_1 + t_frac * uu1_2);
          sample_uu2[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu2_1 + t_frac * uu2_2);
          sample_uu3[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu3_1 + t_frac * uu3_2);
          sample_bb1[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb1_1 + t_frac * bb1_2);
          sample_bb2[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb2_1 + t_frac * bb2_2);
          sample_bb3[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb3_1 + t_frac * bb3_2);
        }
      }

      // Set intrablock interpolated values
      else if (not simulation_block_interp)
      {
        // Extract indices and coefficients
        int b = sample_inds[adaptive_level](m,n,0);
        int k = sample_inds[adaptive_level](m,n,1);
        int j = sample_inds[adaptive_level](m,n,2);
        int i = sample_inds[adaptive_level](m,n,3);
        int t = 0;
        if (slow_light_on)
          t = sample_inds[adaptive_level](m,n,4);
        double f_k = sample_fracs[adaptive_level](m,n,0);
        double f_j = sample_fracs[adaptive_level](m,n,1);
        double f_i = sample_fracs[adaptive_level](m,n,2);

        // Calculate values without temporal interpolation
        if (not (slow_light_on and slow_interp))
        {
          // Perform spatial interpolation
          double rho = InterpolateSimple(grid_prim[t], ind_rho, b, k, j, i, f_k, f_j, f_i);
          double pgas = InterpolateSimple(grid_prim[t], ind_pgas, b, k, j, i, f_k, f_j, f_i);
          double kappa = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa = InterpolateSimple(grid_prim[t], ind_kappa, b, k, j, i, f_k, f_j, f_i);
          double uu1 = InterpolateSimple(grid_prim[t], ind_uu1, b, k, j, i, f_k, f_j, f_i);
          double uu2 = InterpolateSimple(grid_prim[t], ind_uu2, b, k, j, i, f_k, f_j, f_i);
          double uu3 = InterpolateSimple(grid_prim[t], ind_uu3, b, k, j, i, f_k, f_j, f_i);
          double bb1 = InterpolateSimple(grid_bb[t], ind_bb1, b, k, j, i, f_k, f_j, f_i);
          double bb2 = InterpolateSimple(grid_bb[t], ind_bb2, b, k, j, i, f_k, f_j, f_i);
          double bb3 = InterpolateSimple(grid_bb[t], ind_bb3, b, k, j, i, f_k, f_j, f_i);

          // Account for possible invalid values
          if (rho <= 0.0)
            rho = static_cast<double>(grid_prim[t](ind_rho,b,k,j,i));
          if (pgas <= 0.0)
            pgas = static_cast<double>(grid_prim[t](ind_pgas,b,k,j,i));
          if (plasma_model == PlasmaModel::code_kappa and kappa <= 0.0)
            kappa = static_cast<double>(grid_prim[t](ind_kappa,b,k,j,i));

          // Assign values
          sample_rho[adaptive_level](m,n) = static_cast<float>(rho);
          sample_pgas[adaptive_level](m,n) = static_cast<float>(pgas);
          if (plasma_model == PlasmaModel::code_kappa)
            sample_kappa[adaptive_level](m,n) = static_cast<float>(kappa);
          sample_uu1[adaptive_level](m,n) = static_cast<float>(uu1);
          sample_uu2[adaptive_level](m,n) = static_cast<float>(uu2);
          sample_uu3[adaptive_level](m,n) = static_cast<float>(uu3);
          sample_bb1[adaptive_level](m,n) = static_cast<float>(bb1);
          sample_bb2[adaptive_level](m,n) = static_cast<float>(bb2);
          sample_bb3[adaptive_level](m,n) = static_cast<float>(bb3);
        }

        // Calculate values with temporal interpolation
        else
        {
          // Perform spatial interpolation on first slice
          double rho_1 = InterpolateSimple(grid_prim[t], ind_rho, b, k, j, i, f_k, f_j, f_i);
          double pgas_1 = InterpolateSimple(grid_prim[t], ind_pgas, b, k, j, i, f_k, f_j, f_i);
          double kappa_1 = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa_1 = InterpolateSimple(grid_prim[t], ind_kappa, b, k, j, i, f_k, f_j, f_i);
          double uu1_1 = InterpolateSimple(grid_prim[t], ind_uu1, b, k, j, i, f_k, f_j, f_i);
          double uu2_1 = InterpolateSimple(grid_prim[t], ind_uu2, b, k, j, i, f_k, f_j, f_i);
          double uu3_1 = InterpolateSimple(grid_prim[t], ind_uu3, b, k, j, i, f_k, f_j, f_i);
          double bb1_1 = InterpolateSimple(grid_bb[t], ind_bb1, b, k, j, i, f_k, f_j, f_i);
          double bb2_1 = InterpolateSimple(grid_bb[t], ind_bb2, b, k, j, i, f_k, f_j, f_i);
          double bb3_1 = InterpolateSimple(grid_bb[t], ind_bb3, b, k, j, i, f_k, f_j, f_i);

          // Account for possible invalid values
          if (rho_1 <= 0.0)
            rho_1 = static_cast<double>(grid_prim[t](ind_rho,b,k,j,i));
          if (pgas_1 <= 0.0)
            pgas_1 = static_cast<double>(grid_prim[t](ind_pgas,b,k,j,i));
          if (plasma_model == PlasmaModel::code_kappa and kappa_1 <= 0.0)
            kappa_1 = static_cast<double>(grid_prim[t](ind_kappa,b,k,j,i));

          // Perform spatial interpolation on second slice
          double rho_2 = InterpolateSimple(grid_prim[t+1], ind_rho, b, k, j, i, f_k, f_j, f_i);
          double pgas_2 = InterpolateSimple(grid_prim[t+1], ind_pgas, b, k, j, i, f_k, f_j, f_i);
          double kappa_2 = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa_2 = InterpolateSimple(grid_prim[t+1], ind_kappa, b, k, j, i, f_k, f_j, f_i);
          double uu1_2 = InterpolateSimple(grid_prim[t+1], ind_uu1, b, k, j, i, f_k, f_j, f_i);
          double uu2_2 = InterpolateSimple(grid_prim[t+1], ind_uu2, b, k, j, i, f_k, f_j, f_i);
          double uu3_2 = InterpolateSimple(grid_prim[t+1], ind_uu3, b, k, j, i, f_k, f_j, f_i);
          double bb1_2 = InterpolateSimple(grid_bb[t+1], ind_bb1, b, k, j, i, f_k, f_j, f_i);
          double bb2_2 = InterpolateSimple(grid_bb[t+1], ind_bb2, b, k, j, i, f_k, f_j, f_i);
          double bb3_2 = InterpolateSimple(grid_bb[t+1], ind_bb3, b, k, j, i, f_k, f_j, f_i);

          // Account for possible invalid values
          if (rho_2 <= 0.0)
            rho_2 = static_cast<double>(grid_prim[t+1](ind_rho,b,k,j,i));
          if (pgas_2 <= 0.0)
            pgas_2 = static_cast<double>(grid_prim[t+1](ind_pgas,b,k,j,i));
          if (plasma_model == PlasmaModel::code_kappa and kappa_2 <= 0.0)
            kappa_2 = static_cast<double>(grid_prim[t+1](ind_kappa,b,k,j,i));

          // Assign interpolated values
          double t_frac = sample_fracs[adaptive_level](m,n,3);
          sample_rho[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * rho_1 + t_frac * rho_2);
          sample_pgas[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * pgas_1 + t_frac * pgas_2);
          if (plasma_model == PlasmaModel::code_kappa)
            sample_kappa[adaptive_level](m,n) =
                static_cast<float>((1.0 - t_frac) * kappa_1 + t_frac * kappa_2);
          sample_uu1[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu1_1 + t_frac * uu1_2);
          sample_uu2[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu2_1 + t_frac * uu2_2);
          sample_uu3[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu3_1 + t_frac * uu3_2);
          sample_bb1[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb1_1 + t_frac * bb1_2);
          sample_bb2[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb2_1 + t_frac * bb2_2);
          sample_bb3[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb3_1 + t_frac * bb3_2);
        }
      }

      // Set interblock interpolated values
      else
      {
        // Extract index
        int t = 0;
        if (slow_light_on)
          t = sample_inds[adaptive_level](m,n,4);

        // Calculate values without temporal interpolation
        if (not (slow_light_on and slow_interp))
        {
          // Perform spatial interpolation
          double rho = InterpolateAdvanced(grid_prim[t], ind_rho, m, n);
          double pgas = InterpolateAdvanced(grid_prim[t], ind_pgas, m, n);
          double kappa = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa = InterpolateAdvanced(grid_prim[t], ind_kappa, m, n);
          double uu1 = InterpolateAdvanced(grid_prim[t], ind_uu1, m, n);
          double uu2 = InterpolateAdvanced(grid_prim[t], ind_uu2, m, n);
          double uu3 = InterpolateAdvanced(grid_prim[t], ind_uu3, m, n);
          double bb1 = InterpolateAdvanced(grid_bb[t], ind_bb1, m, n);
          double bb2 = InterpolateAdvanced(grid_bb[t], ind_bb2, m, n);
          double bb3 = InterpolateAdvanced(grid_bb[t], ind_bb3, m, n);

          // Account for possible invalid values
          int b = sample_inds[adaptive_level](m,n,0,0);
          int k = sample_inds[adaptive_level](m,n,0,1);
          int j = sample_inds[adaptive_level](m,n,0,2);
          int i = sample_inds[adaptive_level](m,n,0,3);
          if (rho <= 0.0)
            rho = static_cast<double>(grid_prim[t](ind_rho,b,k,j,i));
          if (pgas <= 0.0)
            pgas = static_cast<double>(grid_prim[t](ind_pgas,b,k,j,i));
          if (plasma_model == PlasmaModel::code_kappa and kappa <= 0.0)
            kappa = static_cast<double>(grid_prim[t](ind_kappa,b,k,j,i));

          // Assign values
          sample_rho[adaptive_level](m,n) = static_cast<float>(rho);
          sample_pgas[adaptive_level](m,n) = static_cast<float>(pgas);
          if (plasma_model == PlasmaModel::code_kappa)
            sample_kappa[adaptive_level](m,n) = static_cast<float>(kappa);
          sample_uu1[adaptive_level](m,n) = static_cast<float>(uu1);
          sample_uu2[adaptive_level](m,n) = static_cast<float>(uu2);
          sample_uu3[adaptive_level](m,n) = static_cast<float>(uu3);
          sample_bb1[adaptive_level](m,n) = static_cast<float>(bb1);
          sample_bb2[adaptive_level](m,n) = static_cast<float>(bb2);
          sample_bb3[adaptive_level](m,n) = static_cast<float>(bb3);
        }

        // Calculate values with temporal interpolation
        else
        {
          // Perform spatial interpolation on first slice
          double rho_1 = InterpolateAdvanced(grid_prim[t], ind_rho, m, n);
          double pgas_1 = InterpolateAdvanced(grid_prim[t], ind_pgas, m, n);
          double kappa_1 = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa_1 = InterpolateAdvanced(grid_prim[t], ind_kappa, m, n);
          double uu1_1 = InterpolateAdvanced(grid_prim[t], ind_uu1, m, n);
          double uu2_1 = InterpolateAdvanced(grid_prim[t], ind_uu2, m, n);
          double uu3_1 = InterpolateAdvanced(grid_prim[t], ind_uu3, m, n);
          double bb1_1 = InterpolateAdvanced(grid_bb[t], ind_bb1, m, n);
          double bb2_1 = InterpolateAdvanced(grid_bb[t], ind_bb2, m, n);
          double bb3_1 = InterpolateAdvanced(grid_bb[t], ind_bb3, m, n);

          // Account for possible invalid values
          int b = sample_inds[adaptive_level](m,n,0,0);
          int k = sample_inds[adaptive_level](m,n,0,1);
          int j = sample_inds[adaptive_level](m,n,0,2);
          int i = sample_inds[adaptive_level](m,n,0,3);
          if (rho_1 <= 0.0)
            rho_1 = static_cast<double>(grid_prim[t](ind_rho,b,k,j,i));
          if (pgas_1 <= 0.0)
            pgas_1 = static_cast<double>(grid_prim[t](ind_pgas,b,k,j,i));
          if (plasma_model == PlasmaModel::code_kappa and kappa_1 <= 0.0)
            kappa_1 = static_cast<double>(grid_prim[t](ind_kappa,b,k,j,i));

          // Perform spatial interpolation on second slice
          double rho_2 = InterpolateAdvanced(grid_prim[t+1], ind_rho, m, n);
          double pgas_2 = InterpolateAdvanced(grid_prim[t+1], ind_pgas, m, n);
          double kappa_2 = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa_2 = InterpolateAdvanced(grid_prim[t+1], ind_kappa, m, n);
          double uu1_2 = InterpolateAdvanced(grid_prim[t+1], ind_uu1, m, n);
          double uu2_2 = InterpolateAdvanced(grid_prim[t+1], ind_uu2, m, n);
          double uu3_2 = InterpolateAdvanced(grid_prim[t+1], ind_uu3, m, n);
          double bb1_2 = InterpolateAdvanced(grid_bb[t+1], ind_bb1, m, n);
          double bb2_2 = InterpolateAdvanced(grid_bb[t+1], ind_bb2, m, n);
          double bb3_2 = InterpolateAdvanced(grid_bb[t+1], ind_bb3, m, n);

          // Account for possible invalid values
          if (rho_2 <= 0.0)
            rho_2 = static_cast<double>(grid_prim[t+1](ind_rho,b,k,j,i));
          if (pgas_2 <= 0.0)
            pgas_2 = static_cast<double>(grid_prim[t+1](ind_pgas,b,k,j,i));
          if (plasma_model == PlasmaModel::code_kappa and kappa_2 <= 0.0)
            kappa_2 = static_cast<double>(grid_prim[t+1](ind_kappa,b,k,j,i));

          // Assign interpolated values
          double t_frac = sample_fracs[adaptive_level](m,n,3);
          sample_rho[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * rho_1 + t_frac * rho_2);
          sample_pgas[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * pgas_1 + t_frac * pgas_2);
          if (plasma_model == PlasmaModel::code_kappa)
            sample_kappa[adaptive_level](m,n) =
                static_cast<float>((1.0 - t_frac) * kappa_1 + t_frac * kappa_2);
          sample_uu1[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu1_1 + t_frac * uu1_2);
          sample_uu2[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu2_1 + t_frac * uu2_2);
          sample_uu3[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * uu3_1 + t_frac * uu3_2);
          sample_bb1[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb1_1 + t_frac * bb1_2);
          sample_bb2[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb2_1 + t_frac * bb2_2);
          sample_bb3[adaptive_level](m,n) =
              static_cast<float>((1.0 - t_frac) * bb3_1 + t_frac * bb3_2);
        }
      }
    }
  }

  // Free memory
  if (adaptive_level > 0)
  {
    sample_inds[adaptive_level].Deallocate();
    sample_fracs[adaptive_level].Deallocate();
    sample_nan[adaptive_level].Deallocate();
    sample_fallback[adaptive_level].Deallocate();
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for finding cell indices in a given or nearby block
// Inputs:
//   b: block index
//   k, j, i: cell indices for x^3, x^2, and x^1, possibly 1 beyond valid range
//   k_c, j_c, i_c: cell indices for closest cell to sample point in question
//   x3, x2, x1: coordinates of sample point in question
// Outputs:
//   inds: block, x^3, x^2, and x^1 indices set
// Notes:
//   If requested cell is within block, that cell is used.
//   If requested cell is in an adjacent block, the appropriate cell there is used:
//     If adjacent block is at the same refinement level, the appropriate ghost cell is used.
//     If adjacent block is coarser, the coarse cell containing the appropriate ghost cell is used.
//     If adjacent block is finer, the fine cell closest to the given coordinates is used.
//   If the requested cell is not on the grid, values are copied from the unique cell on the grid
//       closest to the appropriate ghost cell, effectively resulting in constant (rather than
//       linear) extrapolation near the edges of the grid.
//   In the case of simulation_coord == Coordinates::sph_ks, neighboring blocks are understood to
//       cross the periodic boundary in x^3 (phi), but the domain is not stitched together at the
//       poles.
void RadiationIntegrator::FindNearbyInds(int b, int k, int j, int i, int k_c, int j_c, int i_c,
    double x3, double x2, double x1, int inds[4])
{
  // Extract location data
  int n_b = x1f.n2;
  int n_i = x1v.n1;
  int n_j = x2v.n1;
  int n_k = x3v.n1;
  int level = levels(b);
  int location_i = locations(b,0);
  int location_j = locations(b,1);
  int location_k = locations(b,2);
  bool upper_i = i > n_i / 2;
  bool upper_j = j > n_j / 2;
  bool upper_k = k > n_k / 2;
  int i_safe = std::max(std::min(i, n_i - 1), 0);
  int j_safe = std::max(std::min(j, n_j - 1), 0);
  int k_safe = std::max(std::min(k, n_k - 1), 0);

  // Handle simple case on given block
  if (i == i_safe and j == j_safe and k == k_safe)
  {
    inds[0] = b;
    inds[1] = k;
    inds[2] = j;
    inds[3] = i;
    return;
  }

  // Check for grid existing in various directions
  bool x1_off_grid = true;
  bool x2_off_grid = true;
  bool x3_off_grid = true;
  for (int b_alt = 0; b_alt < n_b; b_alt++)
  {
    // Extract data for block
    int level_alt = levels(b_alt);
    int location_i_alt = locations(b_alt,0);
    int location_j_alt = locations(b_alt,1);
    int location_k_alt = locations(b_alt,2);

    // Check x^1-direction
    if (x1_off_grid and i != i_safe) {
      bool same_level_exists = level_alt == level;
      same_level_exists =
          same_level_exists and location_i_alt == (i == -1 ? location_i - 1 : location_i + 1);
      same_level_exists = same_level_exists and location_j_alt == location_j;
      same_level_exists = same_level_exists and location_k_alt == location_k;
      bool coarser_level_exists = level_alt == level - 1;
      coarser_level_exists = coarser_level_exists
          and location_i_alt == (i == -1 ? (location_i - 1) / 2 : (location_i + 1) / 2);
      coarser_level_exists = coarser_level_exists and location_j_alt == location_j / 2;
      coarser_level_exists = coarser_level_exists and location_k_alt == location_k / 2;
      bool finer_level_exists = level_alt == level + 1;
      finer_level_exists = finer_level_exists
          and location_i_alt == (i == -1 ? location_i * 2 - 1 : location_i * 2 + 2);
      finer_level_exists =
          finer_level_exists and location_j_alt == (upper_j ? location_j * 2 + 1 : location_j * 2);
      finer_level_exists =
          finer_level_exists and location_k_alt == (upper_k ? location_k * 2 + 1 : location_k * 2);
      if (same_level_exists or coarser_level_exists or finer_level_exists)
        x1_off_grid = false;
    }

    // Check x^2-direction
    if (x2_off_grid and j != j_safe) {
      bool same_level_exists = level_alt == level;
      same_level_exists = same_level_exists and location_i_alt == location_i;
      same_level_exists =
          same_level_exists and location_j_alt == (j == -1 ? location_j - 1 : location_j + 1);
      same_level_exists = same_level_exists and location_k_alt == location_k;
      bool coarser_level_exists = level_alt == level - 1;
      coarser_level_exists = coarser_level_exists and location_i_alt == location_i / 2;
      coarser_level_exists = coarser_level_exists
          and location_j_alt == (j == -1 ? (location_j - 1) / 2 : (location_j + 1) / 2);
      coarser_level_exists = coarser_level_exists and location_k_alt == location_k / 2;
      bool finer_level_exists = level_alt == level + 1;
      finer_level_exists =
          finer_level_exists and location_i_alt == (upper_i ? location_i * 2 + 1 : location_i * 2);
      finer_level_exists = finer_level_exists
          and location_j_alt == (j == -1 ? location_j * 2 - 1 : location_j * 2 + 2);
      finer_level_exists =
          finer_level_exists and location_k_alt == (upper_k ? location_k * 2 + 1 : location_k * 2);
      if (same_level_exists or coarser_level_exists or finer_level_exists)
        x2_off_grid = false;
    }

    // Check x^3-direction
    if (x3_off_grid and k != k_safe) {
      bool same_level_exists = level_alt == level;
      same_level_exists = same_level_exists and location_i_alt == location_i;
      same_level_exists = same_level_exists and location_j_alt == location_j;
      same_level_exists =
          same_level_exists and location_k_alt == (k == -1 ? location_k - 1 : location_k + 1);
      bool coarser_level_exists = level_alt == level - 1;
      coarser_level_exists = coarser_level_exists and location_i_alt == location_i / 2;
      coarser_level_exists = coarser_level_exists and location_j_alt == location_j / 2;
      coarser_level_exists = coarser_level_exists
          and location_k_alt == (k == -1 ? (location_k - 1) / 2 : (location_k + 1) / 2);
      bool finer_level_exists = level_alt == level + 1;
      finer_level_exists =
          finer_level_exists and location_i_alt == (upper_i ? location_i * 2 + 1 : location_i * 2);
      finer_level_exists =
          finer_level_exists and location_j_alt == (upper_j ? location_j * 2 + 1 : location_j * 2);
      finer_level_exists = finer_level_exists
          and location_k_alt == (k == -1 ? location_k * 2 - 1 : location_k * 2 + 2);
      if (same_level_exists or coarser_level_exists or finer_level_exists)
        x3_off_grid = false;
    }

    // Check x^3-direction across periodic boundary
    if (x3_off_grid and simulation_coord == Coordinates::sph_ks and k == -1 and location_k == 0) {
      bool same_level_exists = level_alt == level;
      same_level_exists = same_level_exists and location_i_alt == location_i;
      same_level_exists = same_level_exists and location_j_alt == location_j;
      same_level_exists = same_level_exists and location_k_alt == n_3_level(level_alt) - 1;
      bool coarser_level_exists = level_alt == level - 1;
      coarser_level_exists = coarser_level_exists and location_i_alt == location_i / 2;
      coarser_level_exists = coarser_level_exists and location_j_alt == location_j / 2;
      coarser_level_exists = coarser_level_exists and location_k_alt == n_3_level(level_alt) - 1;
      bool finer_level_exists = level_alt == level + 1;
      finer_level_exists =
          finer_level_exists and location_i_alt == (upper_i ? location_i * 2 + 1 : location_i * 2);
      finer_level_exists =
          finer_level_exists and location_j_alt == (upper_j ? location_j * 2 + 1 : location_j * 2);
      finer_level_exists = finer_level_exists and location_k_alt == n_3_level(level_alt) - 1;
      if (same_level_exists or coarser_level_exists or finer_level_exists)
        x3_off_grid = false;
    }
    if (x3_off_grid and simulation_coord == Coordinates::sph_ks and k == n_k
        and location_k == n_3_level(level) - 1) {
      bool same_level_exists = level_alt == level;
      same_level_exists = same_level_exists and location_i_alt == location_i;
      same_level_exists = same_level_exists and location_j_alt == location_j;
      same_level_exists = same_level_exists and location_k_alt == 0;
      bool coarser_level_exists = level_alt == level - 1;
      coarser_level_exists = coarser_level_exists and location_i_alt == location_i / 2;
      coarser_level_exists = coarser_level_exists and location_j_alt == location_j / 2;
      coarser_level_exists = coarser_level_exists and location_k_alt == 0;
      bool finer_level_exists = level_alt == level + 1;
      finer_level_exists =
          finer_level_exists and location_i_alt == (upper_i ? location_i * 2 + 1 : location_i * 2);
      finer_level_exists =
          finer_level_exists and location_j_alt == (upper_j ? location_j * 2 + 1 : location_j * 2);
      finer_level_exists = finer_level_exists and location_k_alt == 0;
      if (same_level_exists or coarser_level_exists or finer_level_exists)
        x3_off_grid = false;
    }
  }

  // Account for grid existing in simple cases
  if (i == i_safe)
    x1_off_grid = false;
  if (j == j_safe)
    x2_off_grid = false;
  if (k == k_safe)
    x3_off_grid = false;

  // Adjust sought location to be on grid
  if (x1_off_grid)
    i = i_safe;
  if (x2_off_grid)
    j = j_safe;
  if (x3_off_grid)
    k = k_safe;

  // Find cell at same level
  int level_sought = level;
  int location_i_sought = i == i_safe ? location_i : i == -1 ? location_i - 1 : location_i + 1;
  int location_j_sought = j == j_safe ? location_j : j == -1 ? location_j - 1 : location_j + 1;
  int location_k_sought = k == k_safe ? location_k : k == -1 ? location_k - 1 : location_k + 1;
  if (simulation_coord == Coordinates::sph_ks and k == -1 and location_k == 0)
    location_k_sought = n_3_level(level_sought) - 1;
  if (simulation_coord == Coordinates::sph_ks and k == n_k and location_k == n_3_level(level) - 1)
    location_k_sought = 0;
  int i_sought = i == i_safe ? i : i == -1 ? n_i - 1 : 0;
  int j_sought = j == j_safe ? j : j == -1 ? n_j - 1 : 0;
  int k_sought = k == k_safe ? k : k == -1 ? n_k - 1 : 0;
  for (int b_alt = 0; b_alt < n_b; b_alt++)
    if (levels(b_alt) == level_sought and locations(b_alt,0) == location_i_sought
        and locations(b_alt,1) == location_j_sought and locations(b_alt,2) == location_k_sought)
    {
      inds[0] = b_alt;
      inds[1] = k_sought;
      inds[2] = j_sought;
      inds[3] = i_sought;
      return;
    }

  // Find cell at coarser level
  level_sought = level - 1;
  if (level_sought >= 0)
  {
    location_i_sought =
        i == i_safe ? location_i / 2 : i == -1 ? (location_i - 1) / 2 : (location_i + 1) / 2;
    location_j_sought =
        j == j_safe ? location_j / 2 : j == -1 ? (location_j - 1) / 2 : (location_j + 1) / 2;
    location_k_sought =
        k == k_safe ? location_k / 2 : k == -1 ? (location_k - 1) / 2 : (location_k + 1) / 2;
    if (simulation_coord == Coordinates::sph_ks and k == -1 and location_k == 0)
      location_k_sought = n_3_level(level_sought) - 1;
    if (simulation_coord == Coordinates::sph_ks and k == n_k and location_k == n_3_level(level) - 1)
      location_k_sought = 0;
    i_sought = i == i_safe ? (location_i % 2 * n_i + i) / 2 : i == -1 ? n_i - 1 : 0;
    j_sought = j == j_safe ? (location_j % 2 * n_j + j) / 2 : j == -1 ? n_j - 1 : 0;
    k_sought = k == k_safe ? (location_k % 2 * n_k + k) / 2 : k == -1 ? n_k - 1 : 0;
    for (int b_alt = 0; b_alt < n_b; b_alt++)
      if (levels(b_alt) == level_sought and locations(b_alt,0) == location_i_sought
          and locations(b_alt,1) == location_j_sought and locations(b_alt,2) == location_k_sought)
      {
        inds[0] = b_alt;
        inds[1] = k_sought;
        inds[2] = j_sought;
        inds[3] = i_sought;
        return;
      }
  }

  // Find cell at finer level
  level_sought = level + 1;
  location_i_sought = location_i * 2 + (i == i_safe ? 0 : i == -1 ? -1 : 1) + (upper_i ? 1 : 0);
  location_j_sought = location_j * 2 + (j == j_safe ? 0 : j == -1 ? -1 : 1) + (upper_j ? 1 : 0);
  location_k_sought = location_k * 2 + (k == k_safe ? 0 : k == -1 ? -1 : 1) + (upper_k ? 1 : 0);
  if (simulation_coord == Coordinates::sph_ks and k == -1 and location_k == 0
      and level_sought <= max_level)
    location_k_sought = n_3_level(level_sought) - 1;
  if (simulation_coord == Coordinates::sph_ks and k == n_k and location_k == n_3_level(level) - 1)
    location_k_sought = 0;
  i_sought = i == i_safe ? (upper_i ? (i - n_i / 2) * 2 : i * 2) : i == -1 ? n_i - 2 : 0;
  j_sought = j == j_safe ? (upper_j ? (j - n_j / 2) * 2 : j * 2) : j == -1 ? n_j - 2 : 0;
  k_sought = k == k_safe ? (upper_k ? (k - n_k / 2) * 2 : k * 2) : k == -1 ? n_k - 2 : 0;
  for (int b_alt = 0; b_alt < n_b; b_alt++)
    if (levels(b_alt) == level_sought and locations(b_alt,0) == location_i_sought
        and locations(b_alt,1) == location_j_sought and locations(b_alt,2) == location_k_sought)
    {
      inds[0] = b_alt;
      inds[1] = k_sought;
      inds[2] = j_sought;
      inds[3] = i_sought;
      inds[1] += k < k_c or (k == k_c and x3 > static_cast<double>(x3v(b,k_c))) ? 1 : 0;
      inds[2] += j < j_c or (j == j_c and x2 > static_cast<double>(x2v(b,j_c))) ? 1 : 0;
      inds[3] += i < i_c or (i == i_c and x1 > static_cast<double>(x1v(b,i_c))) ? 1 : 0;
      return;
    }

  // Report grid inconsistency
  throw BlacklightException("Grid interpolation failed.");
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for performing simple interpolation near a given cell
// Inputs:
//   grid_vals: full array of values on grid
//   grid_ind: index of quantity to be interpolated
//   b: block index
//   k, j, i: cell indices
//   f_k, f_j, f_i: interpolation fractions
// Outputs:
//   returned value: interpolated value from grid
double RadiationIntegrator::InterpolateSimple(const Array<float> &grid_vals, int grid_ind, int b,
    int k, int j, int i, double f_k, double f_j, double f_i)
{
  double val_mmm = static_cast<double>(grid_vals(grid_ind,b,k,j,i));
  double val_mmp = static_cast<double>(grid_vals(grid_ind,b,k,j,i+1));
  double val_mpm = static_cast<double>(grid_vals(grid_ind,b,k,j+1,i));
  double val_mpp = static_cast<double>(grid_vals(grid_ind,b,k,j+1,i+1));
  double val_pmm = static_cast<double>(grid_vals(grid_ind,b,k+1,j,i));
  double val_pmp = static_cast<double>(grid_vals(grid_ind,b,k+1,j,i+1));
  double val_ppm = static_cast<double>(grid_vals(grid_ind,b,k+1,j+1,i));
  double val_ppp = static_cast<double>(grid_vals(grid_ind,b,k+1,j+1,i+1));
  double val = (1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
      + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm
      + (1.0 - f_k) * f_j * f_i * val_mpp + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm
      + f_k * (1.0 - f_j) * f_i * val_pmp + f_k * f_j * (1.0 - f_i) * val_ppm
      + f_k * f_j * f_i * val_ppp;
  return val;
}

//--------------------------------------------------------------------------------------------------

// Function for performing advanced interpolation among 8 cells
// Inputs:
//   grid_vals: full array of values on grid
//   grid_ind: index of quantity to be interpolated
//   m: ray index
//   n: index along ray
// Outputs:
//   returned value: interpolated value from grid
// Notes:
//   Assumes sample_inds[adaptive_level] and sample_fracs[adaptive_level] have been set.
double RadiationIntegrator::InterpolateAdvanced(const Array<float> &grid_vals, int grid_ind, int m,
    int n)
{
  double vals[8] = {};
  for (int p = 0; p < 8; p++)
  {
    int b = sample_inds[adaptive_level](m,n,p,0);
    int k = sample_inds[adaptive_level](m,n,p,1);
    int j = sample_inds[adaptive_level](m,n,p,2);
    int i = sample_inds[adaptive_level](m,n,p,3);
    vals[p] = static_cast<double>(grid_vals(grid_ind,b,k,j,i));
  }
  double f_k = sample_fracs[adaptive_level](m,n,0);
  double f_j = sample_fracs[adaptive_level](m,n,1);
  double f_i = sample_fracs[adaptive_level](m,n,2);
  double val = (1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0]
      + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1] + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2]
      + (1.0 - f_k) * f_j * f_i * vals[3] + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4]
      + f_k * (1.0 - f_j) * f_i * vals[5] + f_k * f_j * (1.0 - f_i) * vals[6]
      + f_k * f_j * f_i * vals[7];
  return val;
}
