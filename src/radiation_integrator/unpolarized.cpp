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
//   Works on arrays appropriate for root level or adaptively refined regions.
//   Assumes sample_num (or sample_num_adaptive[adaptive_current_level]), sample_len (or
//       sample_len_adaptive[adaptive_current_level]), j_i (or j_i_adaptive), and alpha_i (or
//       alpha_i_adaptive) have been set.
//   Assumes sample_pos (or sample_pos_adaptive) has been set if image_time == true or
//       image_length == true.
//   Assumes sample_dir (or sample_dir_adaptive) has been set if image_length == true.
//   Assumes cell_values (or cell_values_adaptive) has been set if image_lambda_ave == true or
//       image_emission_ave == true or image_tau_int == true.
//   Allocates and initializes image (or image_adaptive[adaptive_current_level]).
//   Deallocates j_i_adaptive and alpha_i_adaptive.
//   Deallocates cell_values_adaptive if render_num_images <= 0.
void RadiationIntegrator::IntegrateUnpolarizedRadiation()
{
  // Allocate image array
  int num_pix = camera_num_pix;
  if (adaptive_on and adaptive_current_level > 0)
  {
    num_pix = block_counts[adaptive_current_level] * block_num_pix;
    image_adaptive[adaptive_current_level].Allocate(image_num_quantities, num_pix);
    image_adaptive[adaptive_current_level].Zero();
  }
  else if (first_time)
  {
    image.Allocate(image_num_quantities, num_pix);
    image.Zero();
  }
  else
    image.Zero();

  // Alias arrays
  Array<int> sample_num_local = sample_num;
  Array<double> sample_pos_local = sample_pos;
  Array<double> sample_dir_local = sample_dir;
  Array<double> sample_len_local = sample_len;
  Array<double> j_i_local = j_i;
  Array<double> alpha_i_local = alpha_i;
  Array<double> cell_values_local = cell_values;
  Array<double> image_local = image;
  if (adaptive_on and adaptive_current_level > 0)
  {
    sample_num_local = sample_num_adaptive[adaptive_current_level];
    sample_pos_local = sample_pos_adaptive[adaptive_current_level];
    sample_dir_local = sample_dir_adaptive[adaptive_current_level];
    sample_len_local = sample_len_adaptive[adaptive_current_level];
    j_i_local = j_i_adaptive;
    alpha_i_local = alpha_i_adaptive;
    cell_values_local = cell_values_adaptive;
    image_local = image_adaptive[adaptive_current_level];
  }

  // Calculate units
  double x_unit = Physics::gg_msun * mass_msun / (Physics::c * Physics::c);
  double t_unit = x_unit / Physics::c;

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

      // Prepare integrated quantities
      double integrated_lambda = 0.0;
      double integrated_emission = 0.0;

      // Go through samples
      for (int n = 0; n < num_steps; n++)
      {
        // Extract and calculate useful values
        double delta_lambda = sample_len_local(m,n);
        double delta_lambda_cgs = delta_lambda * x_unit / momentum_factor;
        double t_cgs = sample_pos_local(m,n,0) * t_unit;
        double x1 = sample_pos_local(m,n,1);
        double x2 = sample_pos_local(m,n,2);
        double x3 = sample_pos_local(m,n,3);
        double kcov[4];
        kcov[0] = sample_dir_local(m,n,0);
        kcov[1] = sample_dir_local(m,n,1);
        kcov[2] = sample_dir_local(m,n,2);
        kcov[3] = sample_dir_local(m,n,3);
        double j = std::numeric_limits<double>::quiet_NaN();
        if (image_light or image_emission or image_emission_ave)
          j = j_i_local(m,n);
        double alpha = std::numeric_limits<double>::quiet_NaN();
        if (image_light or image_tau or image_tau_int)
          alpha = alpha_i_local(m,n);
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
              image_local(0,m) = exp_neg * (image_local(0,m) + ss * expm1);
            else
              image_local(0,m) = ss;
          }
          else
            image_local(0,m) += j * delta_lambda_cgs;
        }

        // Integrate alternative image quantities
        if (image_time)
          image_local(image_offset_time,m) = std::min(image_local(image_offset_time,m), t_cgs);
        if (image_length)
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
          image_local(image_offset_length,m) += std::sqrt(dl_dlambda_sq) * delta_lambda * x_unit;
        }
        if (image_lambda or image_lambda_ave)
          integrated_lambda += delta_lambda_cgs;
        if (image_emission or image_emission_ave)
          integrated_emission += j * delta_lambda_cgs;
        if (image_tau)
          image_local(image_offset_tau,m) += delta_tau;
        if (image_lambda_ave and not std::isnan(cell_values_local(0,m,n)))
          for (int a = 0; a < CellValues::num_cell_values; a++)
            image_local(image_offset_lambda_ave+a,m) += cell_values_local(a,m,n) * delta_lambda_cgs;
        if (image_emission_ave and not std::isnan(cell_values_local(0,m,n)))
          for (int a = 0; a < CellValues::num_cell_values; a++)
            image_local(image_offset_emission_ave+a,m) +=
                cell_values_local(a,m,n) * j * delta_lambda_cgs;
        if (image_tau_int and not std::isnan(cell_values_local(0,m,n)))
        {
          if (optically_thin)
            for (int a = 0; a < CellValues::num_cell_values; a++)
              image_local(image_offset_tau_int+a,m) = exp_neg
                  * (image_local(image_offset_tau_int+a,m) + cell_values_local(a,m,n) * expm1);
          else
            for (int a = 0; a < CellValues::num_cell_values; a++)
              image_local(image_offset_tau_int+a,m) = cell_values_local(a,m,n);
        }
      }

      // Store integrated quantities
      if (image_lambda)
        image_local(image_offset_lambda,m) = integrated_lambda;
      if (image_emission)
        image_local(image_offset_emission,m) = integrated_emission;

      // Normalize integrated quantities
      if (image_lambda_ave)
        for (int a = 0; a < CellValues::num_cell_values; a++)
          image_local(image_offset_lambda_ave+a,m) /= integrated_lambda;
      if (image_emission_ave)
        for (int a = 0; a < CellValues::num_cell_values; a++)
          image_local(image_offset_emission_ave+a,m) /= integrated_emission;
    }

    // Transform I_nu/nu^3 to I_nu
    if (image_light)
    {
      double image_frequency_cu = image_frequency * image_frequency * image_frequency;
      #pragma omp for schedule(static)
      for (int m = 0; m < num_pix; m++)
        image_local(0,m) *= image_frequency_cu;
    }
  }

  // Free memory
  j_i_adaptive.Deallocate();
  alpha_i_adaptive.Deallocate();
  if (render_num_images <= 0)
    cell_values_adaptive.Deallocate();
  return;
}
