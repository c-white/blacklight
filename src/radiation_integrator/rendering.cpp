// Blacklight radiation integrator - rendering of cell quantities with false colors

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // exp, expm1, sqrt
#include <fstream>    // ifstream
#include <limits>     // numeric_limits
#include <string>     // string

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Physics, enums
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for reading streamline data from files
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Allocates and initializes streamlines and all Arrays contained therein.
void RadiationIntegrator::ReadStreamFiles()
{
  // Allocate array to hold all streamlines
  num_streamlines = 0;
  for (int n_i = 0; n_i < render_num_images; n_i++)
  {
    int num_features = render_num_features[n_i];
    for (int n_f = 0; n_f < num_features; n_f++)
      if (render_types[n_i][n_f] == RenderType::line or render_types[n_i][n_f] == RenderType::tube)
        num_streamlines++;
  }
  streamlines = new Array<double>[num_streamlines];

  // Read streamline data
  for (int n_i = 0, n_s = 0; n_i < render_num_images; n_i++)
  {
    int num_features = render_num_features[n_i];
    for (int n_f = 0; n_f < num_features; n_f++)
    {
      if (not (render_types[n_i][n_f] == RenderType::line
          or render_types[n_i][n_f] == RenderType::tube))
        continue;
      int num_points = 0;
      std::string line;
      std::ifstream stream_data(render_stream_files[n_i][n_f]);
      while (std::getline(stream_data, line))
        num_points++;
      streamlines[n_s].Allocate(num_points, 3);
      stream_data.clear();
      stream_data.seekg(0);
      for (int n_p = 0; n_p < num_points; n_p++)
        stream_data >> streamlines[n_s](n_p,0) >> streamlines[n_s](n_p,1)
            >> streamlines[n_s](n_p,2);
      n_s++;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for rendering false color image
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes sample_num[adaptive_level], sample_pos[adaptive_level], sample_dir[adaptive_level],
//       sample_len[adaptive_level], and cell_values[adaptive_level] have been set.
//   Allocates and initializes render[adaptive_level].
//   Deallocates cell_values[adaptive_level] if adaptive_level > 0.
void RadiationIntegrator::Render()
{
  // Allocate rendering array
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
    render[adaptive_level].Allocate(render_num_images, 3, num_pix);
  render[adaptive_level].Zero();

  // Allocate streamline bookkeeping arrays
  if (first_time and num_streamlines > 0)
  {
    current_stream_distances.Allocate(num_threads, num_streamlines);
    previous_stream_distances.Allocate(num_threads, num_streamlines);
    current_stream_indices.Allocate(num_threads, num_streamlines);
    previous_stream_indices.Allocate(num_threads, num_streamlines);
  }

  // Determine if cell values are needed
  bool cell_values_needed = false;
  for (int n_i = 0; n_i < render_num_images; n_i++)
  {
    int num_features = render_num_features[n_i];
    for (int n_f = 0; n_f < num_features; n_f++)
      if (render_types[n_i][n_f] == RenderType::fill or render_types[n_i][n_f] == RenderType::thresh
          or render_types[n_i][n_f] == RenderType::rise
          or render_types[n_i][n_f] == RenderType::fall)
        cell_values_needed = true;
  }

  // Determine if fills are present
  bool fill_present = false;
  for (int n_i = 0; n_i < render_num_images; n_i++)
  {
    int num_features = render_num_features[n_i];
    for (int n_f = 0; n_f < num_features; n_f++)
      if (render_types[n_i][n_f] == RenderType::fill)
        fill_present = true;
  }

  // Calculate unit
  double x_unit = Physics::gg_msun * mass_msun / (Physics::c * Physics::c);

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    double gcov[4][4];
    double gcon[4][4];

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      // Determine thread number
      int thread = omp_get_thread_num();

      // Extract number of steps
      int num_steps = sample_num[adaptive_level](m);

      // Prepare cell values
      double previous_values[CellValues::num_cell_values];
      for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
        previous_values[n_v] = std::numeric_limits<double>::quiet_NaN();
      double current_values[CellValues::num_cell_values];

      // Prepare distances to streamlines
      for (int n_s = 0; n_s < num_streamlines; n_s++)
      {
        previous_stream_distances(thread,n_s) = std::numeric_limits<double>::infinity();
        previous_stream_indices(thread,n_s) = -1;
      }

      // Go through samples
      for (int n = 0; n < num_steps; n++)
      {
        // Extract geometric values
        double delta_lambda = sample_len[adaptive_level](m,n);
        double x1 = sample_pos[adaptive_level](m,n,1);
        double x2 = sample_pos[adaptive_level](m,n,2);
        double x3 = sample_pos[adaptive_level](m,n,3);
        double kcov[4];
        kcov[0] = sample_dir[adaptive_level](m,n,0);
        kcov[1] = sample_dir[adaptive_level](m,n,1);
        kcov[2] = sample_dir[adaptive_level](m,n,2);
        kcov[3] = sample_dir[adaptive_level](m,n,3);

        // Extract cell values
        if (cell_values_needed)
          for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
            current_values[n_v] = cell_values[adaptive_level](n_v,m,n);

        // Calculate length
        double delta_length = 0.0;
        if (fill_present)
        {
          CovariantGeodesicMetric(x1, x2, x3, gcov);
          ContravariantGeodesicMetric(x1, x2, x3, gcon);
          double temp_a[4] = {};
          for (int a = 1; a < 4; a++)
            for (int mu = 0; mu < 4; mu++)
              temp_a[a] += (gcon[a][mu] - gcon[0][a] * gcon[0][mu] / gcon[0][0]) * kcov[mu];
          double dl_dlambda_sq = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              dl_dlambda_sq += gcov[a][b] * temp_a[a] * temp_a[b];
          delta_length = std::sqrt(dl_dlambda_sq) * delta_lambda * x_unit;
        }

        // Calculate distances to streamlines
        for (int n_s = 0; n_s < num_streamlines; n_s++)
        {
          int num_pieces = streamlines[n_s].n2 - 1;
          double min_distance_sq = std::numeric_limits<double>::infinity();
          int min_index = -1;
          for (int n_p = 0; n_p < num_pieces; n_p++)
          {
            double p1_x1 = streamlines[n_s](n_p,0);
            double p1_x2 = streamlines[n_s](n_p,1);
            double p1_x3 = streamlines[n_s](n_p,2);
            double p2_x1 = streamlines[n_s](n_p+1,0);
            double p2_x2 = streamlines[n_s](n_p+1,1);
            double p2_x3 = streamlines[n_s](n_p+1,2);
            double l12_x1 = p2_x1 - p1_x1;
            double l12_x2 = p2_x2 - p1_x2;
            double l12_x3 = p2_x3 - p1_x3;
            double l10_x1 = x1 - p1_x1;
            double l10_x2 = x2 - p1_x2;
            double l10_x3 = x3 - p1_x3;
            double l20_x1 = x1 - p2_x1;
            double l20_x2 = x2 - p2_x2;
            double l20_x3 = x3 - p2_x3;
            double temp_x1 = l10_x2 * l12_x3 - l10_x3 * l12_x2;
            double temp_x2 = l10_x3 * l12_x1 - l10_x1 * l12_x3;
            double temp_x3 = l10_x1 * l12_x2 - l10_x2 * l12_x1;
            double distance_sq = (temp_x1 * temp_x1 + temp_x2 * temp_x2 + temp_x3 * temp_x3)
                / (l12_x1 * l12_x1 + l12_x2 * l12_x2 + l12_x3 * l12_x3);
            if (l10_x1 * l12_x1 + l10_x2 * l12_x2 + l10_x3 * l12_x3 < 0.0)
              distance_sq = l10_x1 * l10_x1 + l10_x2 * l10_x2 + l10_x3 * l10_x3;
            else if (l20_x1 * l12_x1 + l20_x2 * l12_x2 + l20_x3 * l12_x3 > 0.0)
              distance_sq = l20_x1 * l20_x1 + l20_x2 * l20_x2 + l20_x3 * l20_x3;
            if (distance_sq < min_distance_sq)
            {
              min_distance_sq = distance_sq;
              min_index = n_p;
            }
          }
          current_stream_distances(thread,n_s) = min_distance_sq;
          current_stream_indices(thread,n_s) = min_index;
        }

        // Go through rendering images
        for (int n_i = 0, n_s = 0; n_i < render_num_images; n_i++)
        {
          // Extract number of features
          int num_features = render_num_features[n_i];

          // Go through features
          for (int n_f = 0; n_f < num_features; n_f++)
          {
            // Process cell-value fills
            if (render_types[n_i][n_f] == RenderType::fill)
            {
              // Extract relevant quantity
              int n_v = render_quantities[n_i][n_f];
              double current_value = current_values[n_v];

              // Calculate effect of passing through filling region
              if (render_types[n_i][n_f] == RenderType::fill
                  and current_value >= render_min_vals[n_i][n_f]
                  and current_value <= render_max_vals[n_i][n_f])
              {
                double delta_tau = delta_length / render_tau_scales[n_i][n_f];
                bool optically_thin = delta_tau <= delta_tau_max;
                if (optically_thin)
                {
                  double exp_neg = std::exp(-delta_tau);
                  double expm1 = std::expm1(delta_tau);
                  render[adaptive_level](n_i,0,m) =
                      exp_neg * (render[adaptive_level](n_i,0,m) + render_x_vals[n_i][n_f] * expm1);
                  render[adaptive_level](n_i,1,m) =
                      exp_neg * (render[adaptive_level](n_i,1,m) + render_y_vals[n_i][n_f] * expm1);
                  render[adaptive_level](n_i,2,m) =
                      exp_neg * (render[adaptive_level](n_i,2,m) + render_z_vals[n_i][n_f] * expm1);
                }
                else
                {
                  render[adaptive_level](n_i,0,m) = render_x_vals[n_i][n_f];
                  render[adaptive_level](n_i,1,m) = render_y_vals[n_i][n_f];
                  render[adaptive_level](n_i,2,m) = render_z_vals[n_i][n_f];
                }
              }
            }

            // Process cell-value thresholds
            else if (render_types[n_i][n_f] == RenderType::thresh
                or render_types[n_i][n_f] == RenderType::rise
                or render_types[n_i][n_f] == RenderType::fall)
            {
              // Extract relevant quantity
              int n_v = render_quantities[n_i][n_f];
              double previous_value = previous_values[n_v];
              double current_value = current_values[n_v];

              // Determine if threshold has been crossed
              bool threshold_crossed = false;
              bool rise_search = render_types[n_i][n_f] == RenderType::thresh
                  or render_types[n_i][n_f] == RenderType::rise;
              if (rise_search and previous_value < render_thresh_vals[n_i][n_f]
                  and current_value >= render_thresh_vals[n_i][n_f])
                threshold_crossed = true;
              bool fall_search = render_types[n_i][n_f] == RenderType::thresh
                  or render_types[n_i][n_f] == RenderType::fall;
              if (fall_search and previous_value > render_thresh_vals[n_i][n_f]
                  and current_value <= render_thresh_vals[n_i][n_f])
                threshold_crossed = true;

              // Calculate effect of crossing threshold
              if (threshold_crossed)
              {
                double opacity = render_opacities[n_i][n_f];
                render[adaptive_level](n_i,0,m) = (1.0 - opacity) * render[adaptive_level](n_i,0,m)
                    + opacity * render_x_vals[n_i][n_f];
                render[adaptive_level](n_i,1,m) = (1.0 - opacity) * render[adaptive_level](n_i,1,m)
                    + opacity * render_y_vals[n_i][n_f];
                render[adaptive_level](n_i,2,m) = (1.0 - opacity) * render[adaptive_level](n_i,2,m)
                    + opacity * render_z_vals[n_i][n_f];
              }
            }

            // Process streamlines
            else if (render_types[n_i][n_f] == RenderType::line
                or render_types[n_i][n_f] == RenderType::tube)
            {
              double current_min_distance_sq = current_stream_distances(thread,n_s);
              double previous_min_distance_sq = previous_stream_distances(thread,n_s);
              double r_sq = render_r_vals[n_i][n_f] * render_r_vals[n_i][n_f];
              if (previous_min_distance_sq < r_sq and current_min_distance_sq >= r_sq)
              {
                // Account for diffuse lighting of tubes
                double reflectance = 1.0;
                if (render_types[n_i][n_f] == RenderType::tube)
                {
                  double previous_x1 = sample_pos[adaptive_level](m,n-1,1);
                  double previous_x2 = sample_pos[adaptive_level](m,n-1,2);
                  double previous_x3 = sample_pos[adaptive_level](m,n-1,3);
                  double ray_dx1 = x1 - previous_x1;
                  double ray_dx2 = x2 - previous_x2;
                  double ray_dx3 = x3 - previous_x3;
                  int previous_min_index = previous_stream_indices(thread,n_s);
                  double tube_dx1 = streamlines[n_s](previous_min_index+1,0)
                      - streamlines[n_s](previous_min_index,0);
                  double tube_dx2 = streamlines[n_s](previous_min_index+1,1)
                      - streamlines[n_s](previous_min_index,1);
                  double tube_dx3 = streamlines[n_s](previous_min_index+1,2)
                      - streamlines[n_s](previous_min_index,2);
                  double temp =
                      std::sqrt(tube_dx1 * tube_dx1 + tube_dx2 * tube_dx2 + tube_dx3 * tube_dx3);
                  tube_dx1 /= temp;
                  tube_dx2 /= temp;
                  tube_dx3 /= temp;
                  temp = ray_dx1 * tube_dx1 + ray_dx2 * tube_dx2 + ray_dx3 * tube_dx3;
                  double ray_vec_1 = ray_dx1 - temp * tube_dx1;
                  double ray_vec_2 = ray_dx2 - temp * tube_dx2;
                  double ray_vec_3 = ray_dx3 - temp * tube_dx3;
                  double alt_dx1 = previous_x1 - streamlines[n_s](previous_min_index,0);
                  double alt_dx2 = previous_x2 - streamlines[n_s](previous_min_index,1);
                  double alt_dx3 = previous_x3 - streamlines[n_s](previous_min_index,2);
                  temp = alt_dx1 * tube_dx1 + alt_dx2 * tube_dx2 + alt_dx3 * tube_dx3;
                  double alt_vec_1 = alt_dx1 - temp * tube_dx1;
                  double alt_vec_2 = alt_dx2 - temp * tube_dx2;
                  double alt_vec_3 = alt_dx3 - temp * tube_dx3;
                  double ray_vec_sq =
                      ray_vec_1 * ray_vec_1 + ray_vec_2 * ray_vec_2 + ray_vec_3 * ray_vec_3;
                  double alt_vec_sq =
                      alt_vec_1 * alt_vec_1 + alt_vec_2 * alt_vec_2 + alt_vec_3 * alt_vec_3;
                  temp = (ray_vec_1 * alt_vec_1 + ray_vec_2 * alt_vec_2 + ray_vec_3 * alt_vec_3)
                      / ray_vec_sq;
                  double frac = std::sqrt(temp * temp + (r_sq - alt_vec_sq) / ray_vec_sq) - temp;
                  double norm_1 = ray_vec_1 * frac + alt_vec_1;
                  double norm_2 = ray_vec_2 * frac + alt_vec_2;
                  double norm_3 = ray_vec_3 * frac + alt_vec_3;
                  temp = std::sqrt(norm_1 * norm_1 + norm_2 * norm_2 + norm_3 * norm_3);
                  norm_1 /= temp;
                  norm_2 /= temp;
                  norm_3 /= temp;
                  double light_1 = render_lx_vals[n_i][n_f];
                  double light_2 = render_ly_vals[n_i][n_f];
                  double light_3 = render_lz_vals[n_i][n_f];
                  temp = std::sqrt(light_1 * light_1 + light_2 * light_2 + light_3 * light_3);
                  light_1 /= temp;
                  light_2 /= temp;
                  light_3 /= temp;
                  double ambient = render_ambient[n_i][n_f];
                  double diffuse = render_diffuse[n_i][n_f];
                  reflectance =
                      ambient + diffuse * (light_1 * norm_1 + light_2 * norm_2 + light_3 * norm_3);
                }

                // Calculate effect of passing through streamline
                double opacity = render_opacities[n_i][n_f];
                render[adaptive_level](n_i,0,m) = (1.0 - opacity) * render[adaptive_level](n_i,0,m)
                    + opacity * reflectance * render_x_vals[n_i][n_f];
                render[adaptive_level](n_i,1,m) = (1.0 - opacity) * render[adaptive_level](n_i,1,m)
                    + opacity * reflectance * render_y_vals[n_i][n_f];
                render[adaptive_level](n_i,2,m) = (1.0 - opacity) * render[adaptive_level](n_i,2,m)
                    + opacity * reflectance * render_z_vals[n_i][n_f];
              }
              n_s++;
            }
          }
        }

        // Store current values
        for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
          previous_values[n_v] = current_values[n_v];
        for (int n_s = 0; n_s < num_streamlines; n_s++)
        {
          previous_stream_distances(thread,n_s) = current_stream_distances(thread,n_s);
          previous_stream_indices(thread,n_s) = current_stream_indices(thread,n_s);
        }
      }
    }
  }

  // Free memory
  if (adaptive_level > 0)
    cell_values[adaptive_level].Deallocate();
  return;
}
