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
//   Works on arrays appropriate for root level or adaptively refined regions.
//   Assumes sample_num (or sample_num_adaptive[adaptive_current_level]), sample_pos (or
//       sample_pos_adaptive[adaptive_current_level]), sample_dir (or
//       sample_dir_adaptive[adaptive_current_level]), sample_len (or
//       sample_len_adaptive[adaptive_current_level]), and cell_values (or cell_values_adaptive)
//       have been set.
//   Allocates and initializes render (or render_adaptive[adaptive_current_level]).
//   Deallocates cell_values_adaptive.
void RadiationIntegrator::Render()
{
  // Allocate rendering array
  int num_pix = camera_num_pix;
  if (adaptive_on and adaptive_current_level > 0)
  {
    num_pix = block_counts[adaptive_current_level] * block_num_pix;
    render_adaptive[adaptive_current_level].Allocate(render_num_images, 3, num_pix);
    render_adaptive[adaptive_current_level].Zero();
  }
  else if (first_time)
  {
    render.Allocate(render_num_images, 3, num_pix);
    render.Zero();
  }
  else
    render.Zero();

  // Alias arrays
  Array<int> sample_num_local = sample_num;
  Array<double> sample_pos_local = sample_pos;
  Array<double> sample_dir_local = sample_dir;
  Array<double> sample_len_local = sample_len;
  Array<double> cell_values_local = cell_values;
  Array<double> render_local = render;
  if (adaptive_on and adaptive_current_level > 0)
  {
    sample_num_local = sample_num_adaptive[adaptive_current_level];
    sample_pos_local = sample_pos_adaptive[adaptive_current_level];
    sample_dir_local = sample_dir_adaptive[adaptive_current_level];
    sample_len_local = sample_len_adaptive[adaptive_current_level];
    cell_values_local = cell_values_adaptive;
    render_local = render_adaptive[adaptive_current_level];
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
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      // Extract number of steps
      int num_steps = sample_num_local(m);

      // Prepare cell values
      double previous_values[CellValues::num_cell_values];
      for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
        previous_values[n_v] = std::numeric_limits<double>::quiet_NaN();
      double current_values[CellValues::num_cell_values];

      // Go through samples
      for (int n = 0; n < num_steps; n++)
      {
        // Extract useful values
        double delta_lambda = sample_len_local(m,n);
        double x1 = sample_pos_local(m,n,1);
        double x2 = sample_pos_local(m,n,2);
        double x3 = sample_pos_local(m,n,3);
        double kcov[4];
        kcov[0] = sample_dir_local(m,n,0);
        kcov[1] = sample_dir_local(m,n,1);
        kcov[2] = sample_dir_local(m,n,2);
        kcov[3] = sample_dir_local(m,n,3);
        for (int n_v = 0; n_v < CellValues::num_cell_values; n_v++)
          current_values[n_v] = cell_values_local(n_v,m,n);

        // Calculate length
        double delta_length = 0.0;
        if (fill_present)
        {
          CovariantGeodesicMetric(x1, x2, x3, gcov);
          ContravariantGeodesicMetric(x1, x2, x3, gcon);
          double dl_dlambda_sq = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              for (int mu = 0; mu < 4; mu++)
                for (int nu = 0; nu < 4; nu++)
                  dl_dlambda_sq += gcov(a,b) * (gcon(a,mu) - gcon(0,a) * gcon(0,mu) / gcon(0,0))
                      * (gcon(b,nu) - gcon(0,b) * gcon(0,nu) / gcon(0,0)) * kcov[mu] * kcov[nu];
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

            // Calculate effect of crossing threshold
            if ((render_types[n_i][n_f] == RenderType::rise
                and previous_value < render_thresh_vals[n_i][n_f]
                and current_value >= render_thresh_vals[n_i][n_f])
                or (render_types[n_i][n_f] == RenderType::fall
                and previous_value > render_thresh_vals[n_i][n_f]
                and current_value <= render_thresh_vals[n_i][n_f]))
            {
              double opacity = render_opacities[n_i][n_f];
              render_local(n_i,0,m) =
                  (1.0 - opacity) * render_local(n_i,0,m) + opacity * render_x_vals[n_i][n_f];
              render_local(n_i,1,m) =
                  (1.0 - opacity) * render_local(n_i,1,m) + opacity * render_y_vals[n_i][n_f];
              render_local(n_i,2,m) =
                  (1.0 - opacity) * render_local(n_i,2,m) + opacity * render_z_vals[n_i][n_f];
            }

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
                render_local(n_i,0,m) =
                    exp_neg * (render_local(n_i,0,m) + render_x_vals[n_i][n_f] * expm1);
                render_local(n_i,1,m) =
                    exp_neg * (render_local(n_i,1,m) + render_y_vals[n_i][n_f] * expm1);
                render_local(n_i,2,m) =
                    exp_neg * (render_local(n_i,2,m) + render_z_vals[n_i][n_f] * expm1);
              }
              else
              {
                render_local(n_i,0,m) = render_x_vals[n_i][n_f];
                render_local(n_i,1,m) = render_y_vals[n_i][n_f];
                render_local(n_i,2,m) = render_z_vals[n_i][n_f];
              }
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
  cell_values_adaptive.Deallocate();
  return;
}
