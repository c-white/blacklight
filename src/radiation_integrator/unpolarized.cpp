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
// Output: (none)
// Notes:
//   Assumes sample_num, sample_len, j_i, and alpha_i have been set.
//   Allocates and initializes image.
void RadiationIntegrator::IntegrateUnpolarizedRadiation()
{
  // Allocate image array
  image.Allocate(image_resolution, image_resolution);
  image.Zero();

  // Calculate unit
  double x_unit = physics::gg_msun * mass_msun / (physics::c * physics::c);

  // Work in parallel
  #pragma omp parallel
  {
    // Go through pixels and samples
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        int num_steps = sample_num(m,l);
        for (int n = 0; n < num_steps; n++)
        {
          double delta_lambda = sample_len(m,l,n);
          double delta_lambda_cgs = delta_lambda * x_unit / momentum_factor;
          if (alpha_i(m,l,n) > 0.0)
          {
            double delta_tau_nu = alpha_i(m,l,n) * delta_lambda_cgs;
            double ss = j_i(m,l,n) / alpha_i(m,l,n);
            if (delta_tau_nu <= delta_tau_max)
              image(m,l) = std::exp(-delta_tau_nu) * (image(m,l) + ss * std::expm1(delta_tau_nu));
            else
              image(m,l) = ss;
          }
          else
            image(m,l) += j_i(m,l,n) * delta_lambda_cgs;
        }
      }

    // Transform I_nu/nu^3 to I_nu
    double image_frequency_cu = image_frequency * image_frequency * image_frequency;
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
        image(m,l) *= image_frequency_cu;
  }
  return;
}
