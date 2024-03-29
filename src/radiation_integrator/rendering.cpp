// Blacklight radiation integrator - rendering of cell quantities with false colors

// C++ headers
#include <cmath>   // exp, expm1
#include <limits>  // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Physics, enums
#include "../utils/array.hpp"        // Array

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
      // Extract number of steps
      int num_steps = sample_num[adaptive_level](m);

      // Prepare cell values
      double previous_values[CellValues::num_cell_values];
      for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
        previous_values[n_v] = std::numeric_limits<double>::quiet_NaN();
      double current_values[CellValues::num_cell_values];

      // Go through samples
      for (int n = 0; n < num_steps; n++)
      {
        // Extract useful values
        double delta_lambda = sample_len[adaptive_level](m,n);
        double x1 = sample_pos[adaptive_level](m,n,1);
        double x2 = sample_pos[adaptive_level](m,n,2);
        double x3 = sample_pos[adaptive_level](m,n,3);
        double kcov[4];
        kcov[0] = sample_dir[adaptive_level](m,n,0);
        kcov[1] = sample_dir[adaptive_level](m,n,1);
        kcov[2] = sample_dir[adaptive_level](m,n,2);
        kcov[3] = sample_dir[adaptive_level](m,n,3);
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

        // Go through rendering images
        for (int n_i = 0; n_i < render_num_images; n_i++)
        {
          // Extract number of features
          int num_features = render_num_features[n_i];

          // Go through features
          for (int n_f = 0; n_f < num_features; n_f++)
          {
            // Extract relevant quantity
            int n_v = render_quantities[n_i][n_f];
            double previous_value = previous_values[n_v];
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
        }

        // Store current values
        for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
          previous_values[n_v] = current_values[n_v];
      }
    }
  }

  // Free memory
  if (adaptive_level > 0)
    cell_values[adaptive_level].Deallocate();
  return;
}
