// Blacklight radiation integrator - unpolarized radiation integration

// C++ headers
#include <algorithm>  // min
#include <cmath>      // exp, expm1
#include <limits>     // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Physics, enums
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for integrating unpolarized radiative transfer equation
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes sample_num[adaptive_level], sample_len[adaptive_level], j_i[adaptive_level],
//       alpha_i[adaptive_level], and momentum_factors[adaptive_level] have been set.
//   Assumes sample_pos[adaptive_level] has been set if image_time == true or image_length == true.
//   Assumes sample_dir[adaptive_level] has been set if image_length == true.
//   Assumes cell_values[adaptive_level] has been set if image_lambda_ave == true or
//       image_emission_ave == true or image_tau_int == true.
//   Allocates and initializes image[adaptive_level].
//   Deallocates j_i[adaptive_level] and alpha_i[adaptive_level] if adaptive_level > 0.
//   Deallocates cell_values[adaptive_level] if render_num_images <= 0 and adaptive_level > 0.
void RadiationIntegrator::IntegrateUnpolarizedRadiation()
{
  // Allocate image array
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
    image[adaptive_level].Allocate(image_num_quantities, num_pix);
  image[adaptive_level].Zero();

  // Calculate units
  double x_unit = Physics::gg_msun * mass_msun / (Physics::c * Physics::c);
  double t_unit = x_unit / Physics::c;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    double gcov[4][4];
    double gcon[4][4];

    // Go through frequencies and pixels
    #pragma omp for schedule(static) collapse(2)
    for (int l = 0; l < image_num_frequencies; l++)
      for (int m = 0; m < num_pix; m++)
      {
        // Extract number of steps
        int num_steps = sample_num[adaptive_level](m);

        // Prepare integrated quantities
        double integrated_lambda = 0.0;
        double integrated_emission = 0.0;

        // Go through samples
        for (int n = 0; n < num_steps; n++)
        {
          // Extract and calculate useful values
          double delta_lambda = sample_len[adaptive_level](m,n);
          double delta_lambda_cgs =
              delta_lambda * x_unit / (image_frequencies(l) * momentum_factors[adaptive_level](m));
          double t_cgs = sample_pos[adaptive_level](m,n,0) * t_unit;
          double x1 = sample_pos[adaptive_level](m,n,1);
          double x2 = sample_pos[adaptive_level](m,n,2);
          double x3 = sample_pos[adaptive_level](m,n,3);
          double kcov[4];
          kcov[0] = sample_dir[adaptive_level](m,n,0);
          kcov[1] = sample_dir[adaptive_level](m,n,1);
          kcov[2] = sample_dir[adaptive_level](m,n,2);
          kcov[3] = sample_dir[adaptive_level](m,n,3);
          double j = std::numeric_limits<double>::quiet_NaN();
          if (image_light or image_emission or image_emission_ave)
            j = j_i[adaptive_level](l,m,n);
          double alpha = std::numeric_limits<double>::quiet_NaN();
          if (image_light or image_tau or image_tau_int)
            alpha = alpha_i[adaptive_level](l,m,n);
          double ss = j / alpha;
          double delta_tau = alpha * delta_lambda_cgs;
          double exp_neg = std::exp(-delta_tau);
          double expm1 = std::expm1(delta_tau);
          bool optically_thin = delta_tau <= delta_tau_max;

          // Integrate light
          if (image_light)
          {
            if (alpha > 0.0)
            {
              if (optically_thin)
                image[adaptive_level](l,m) = exp_neg * (image[adaptive_level](l,m) + ss * expm1);
              else
                image[adaptive_level](l,m) = ss;
            }
            else
              image[adaptive_level](l,m) += j * delta_lambda_cgs;
          }

          // Integrate alternative image quantities
          if (image_time and l == 0)
            image[adaptive_level](image_offset_time,m) =
                std::min(image[adaptive_level](image_offset_time,m), t_cgs);
          if (image_length and l == 0)
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
            image[adaptive_level](image_offset_length,m) +=
                std::sqrt(dl_dlambda_sq) * delta_lambda * x_unit;
          }
          if (image_lambda or image_lambda_ave)
            integrated_lambda += delta_lambda_cgs;
          if (image_emission or image_emission_ave)
            integrated_emission += j * delta_lambda_cgs;
          if (image_tau)
            image[adaptive_level](image_offset_tau+l,m) += delta_tau;
          if (image_lambda_ave and not std::isnan(cell_values[adaptive_level](0,m,n)))
            for (int a = 0; a < CellValues::num_cell_values; a++)
            {
              int index = image_offset_lambda_ave + l * CellValues::num_cell_values + a;
              image[adaptive_level](index,m) +=
                  cell_values[adaptive_level](a,m,n) * delta_lambda_cgs;
            }
          if (image_emission_ave and not std::isnan(cell_values[adaptive_level](0,m,n)))
            for (int a = 0; a < CellValues::num_cell_values; a++)
            {
              int index = image_offset_emission_ave + l * CellValues::num_cell_values + a;
              image[adaptive_level](index,m) +=
                  cell_values[adaptive_level](a,m,n) * j * delta_lambda_cgs;
            }
          if (image_tau_int and not std::isnan(cell_values[adaptive_level](0,m,n)))
          {
            if (optically_thin)
              for (int a = 0; a < CellValues::num_cell_values; a++)
              {
                int index = image_offset_tau_int + l * CellValues::num_cell_values + a;
                image[adaptive_level](index,m) = exp_neg
                    * (image[adaptive_level](index,m) + cell_values[adaptive_level](a,m,n) * expm1);
              }
            else
              for (int a = 0; a < CellValues::num_cell_values; a++)
              {
                int index = image_offset_tau_int + l * CellValues::num_cell_values + a;
                image[adaptive_level](index,m) = cell_values[adaptive_level](a,m,n);
              }
          }
        }

        // Store integrated quantities
        if (image_lambda)
          image[adaptive_level](image_offset_lambda+l,m) = integrated_lambda;
        if (image_emission)
          image[adaptive_level](image_offset_emission+l,m) = integrated_emission;

        // Normalize integrated quantities
        if (image_lambda_ave)
          for (int a = 0; a < CellValues::num_cell_values; a++)
          {
            int index = image_offset_lambda_ave + l * CellValues::num_cell_values + a;
            image[adaptive_level](index,m) /= integrated_lambda;
          }
        if (image_emission_ave)
          for (int a = 0; a < CellValues::num_cell_values; a++)
          {
            int index = image_offset_emission_ave + l * CellValues::num_cell_values + a;
            image[adaptive_level](index,m) /= integrated_emission;
          }
      }

    // Transform I_nu/nu^3 to I_nu
    if (image_light)
    {
      #pragma omp for schedule(static) collapse(2)
      for (int l = 0; l < image_num_frequencies; l++)
        for (int m = 0; m < num_pix; m++)
        {
          double nu_cu = image_frequencies(l) * image_frequencies(l) * image_frequencies(l);
          image[adaptive_level](l,m) *= nu_cu;
        }
    }
  }

  // Free memory
  if (adaptive_level > 0)
  {
    j_i[adaptive_level].Deallocate();
    alpha_i[adaptive_level].Deallocate();
    if (render_num_images <= 0)
      cell_values[adaptive_level].Deallocate();
  }
  return;
}
