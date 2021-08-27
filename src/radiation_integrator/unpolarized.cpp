// Blacklight radiation integrator - unpolarized radiation integration

// C++ headers
#include <cmath>  // exp, expm1

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // physics
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
//   Allocates and initializes image (or image_adaptive[adaptive_current_level]).
//   Deallocates j_i_adaptive and alpha_i_adaptive.
void RadiationIntegrator::IntegrateUnpolarizedRadiation()
{
  // Allocate image array
  int num_pix = camera_num_pix;
  if (adaptive_on and adaptive_current_level > 0)
  {
    num_pix = block_counts[adaptive_current_level] * block_num_pix;
    image_adaptive[adaptive_current_level].Allocate(num_pix);
    image_adaptive[adaptive_current_level].Zero();
  }
  else if (first_time)
  {
    image.Allocate(num_pix);
    image.Zero();
  }

  // Alias arrays
  Array<int> sample_num_local = sample_num;
  Array<double> sample_len_local = sample_len;
  Array<double> j_i_local = j_i;
  Array<double> alpha_i_local = alpha_i;
  Array<double> image_local = image;
  if (adaptive_on and adaptive_current_level > 0)
  {
    sample_num_local = sample_num_adaptive[adaptive_current_level];
    sample_len_local = sample_len_adaptive[adaptive_current_level];
    j_i_local = j_i_adaptive;
    alpha_i_local = alpha_i_adaptive;
    image_local = image_adaptive[adaptive_current_level];
  }

  // Calculate unit
  double x_unit = physics::gg_msun * mass_msun / (physics::c * physics::c);

  // Work in parallel
  #pragma omp parallel
  {
    // Go through pixels and samples
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      int num_steps = sample_num_local(m);
      for (int n = 0; n < num_steps; n++)
      {
        double delta_lambda = sample_len_local(m,n);
        double delta_lambda_cgs = delta_lambda * x_unit / momentum_factor;
        if (alpha_i_local(m,n) > 0.0)
        {
          double delta_tau_nu = alpha_i_local(m,n) * delta_lambda_cgs;
          double ss = j_i_local(m,n) / alpha_i_local(m,n);
          if (delta_tau_nu <= delta_tau_max)
            image_local(m) =
                std::exp(-delta_tau_nu) * (image_local(m) + ss * std::expm1(delta_tau_nu));
          else
            image_local(m) = ss;
        }
        else
          image_local(m) += j_i_local(m,n) * delta_lambda_cgs;
      }
    }

    // Transform I_nu/nu^3 to I_nu
    double image_frequency_cu = image_frequency * image_frequency * image_frequency;
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
      image_local(m) *= image_frequency_cu;
  }

  // Free memory
  j_i_adaptive.Deallocate();
  alpha_i_adaptive.Deallocate();
  return;
}
