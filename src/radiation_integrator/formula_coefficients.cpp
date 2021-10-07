// Blacklight radiation integrator - formula radiative transfer coefficients

// C++ headers
#include <cmath>   // atan, atan2, cos, exp, pow, sin, sqrt
#include <limits>  // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on formula
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Works on arrays appropriate for root level or adaptively refined regions.
//   Assumes sample_flags (or sample_flags_adaptive[adaptive_current_level]), sample_num (or
//       sample_num_adaptive[adaptive_current_level]), sample_pos (or
//       sample_pos_adaptive[adaptive_current_level]), and sample_dir (or
//       sample_dir_adaptive[adaptive_current_level]) have been set.
//   Allocates and initializes j_i (or j_i_adaptive) and alpha_i (or alpha_i_adaptive).
//   References code comparison paper 2020 ApJ 897 148 (C).
void RadiationIntegrator::CalculateFormulaCoefficients()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_on and adaptive_current_level > 0)
  {
    num_pix = block_counts[adaptive_current_level] * block_num_pix;
    int geodesic_num_steps_local = geodesic_num_steps_adaptive[adaptive_current_level];
    j_i_adaptive.Allocate(num_pix, geodesic_num_steps_local);
    alpha_i_adaptive.Allocate(num_pix, geodesic_num_steps_local);
    j_i_adaptive.Zero();
    alpha_i_adaptive.Zero();
  }
  else if (first_time)
  {
    j_i.Allocate(num_pix, geodesic_num_steps);
    alpha_i.Allocate(num_pix, geodesic_num_steps);
    j_i.Zero();
    alpha_i.Zero();
  }
  else
  {
    j_i.Zero();
    alpha_i.Zero();
  }

  // Alias arrays
  Array<bool> sample_flags_local = sample_flags;
  Array<int> sample_num_local = sample_num;
  Array<double> sample_pos_local = sample_pos;
  Array<double> sample_dir_local = sample_dir;
  Array<double> j_i_local = j_i;
  Array<double> alpha_i_local = alpha_i;
  if (adaptive_on and adaptive_current_level > 0)
  {
    sample_flags_local = sample_flags_adaptive[adaptive_current_level];
    sample_num_local = sample_num_adaptive[adaptive_current_level];
    sample_pos_local = sample_pos_adaptive[adaptive_current_level];
    sample_dir_local = sample_dir_adaptive[adaptive_current_level];
    j_i_local = j_i_adaptive;
    alpha_i_local = alpha_i_adaptive;
  }

  // Go through rays in parallel
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < num_pix; m++)
  {
    // Check number of steps
    int num_steps = sample_num_local(m);
    if (num_steps <= 0)
      continue;

    // Set pixel to NaN if ray has problem
    if (fallback_nan and sample_flags_local(m))
    {
      for (int n = 0; n < num_steps; n++)
      {
        j_i_local(m,n) = std::numeric_limits<double>::quiet_NaN();
        alpha_i_local(m,n) = std::numeric_limits<double>::quiet_NaN();
      }
      continue;
    }

    // Go through samples
    for (int n = 0; n < num_steps; n++)
    {
      // Extract geodesic position and momentum
      double x = sample_pos_local(m,n,1);
      double y = sample_pos_local(m,n,2);
      double z = sample_pos_local(m,n,3);
      double k_0 = sample_dir_local(m,n,0);
      double k_1 = sample_dir_local(m,n,1);
      double k_2 = sample_dir_local(m,n,2);
      double k_3 = sample_dir_local(m,n,3);

      // Calculate coordinates and skip coupling if outside camera radius
      double r = RadialGeodesicCoordinate(x, y, z);
      if (r > camera_r)
        continue;
      double rr = std::sqrt(r * r - z * z);
      double cth = z / r;
      double sth = std::sqrt(1.0 - cth * cth);
      double ph = std::atan2(y, x) - std::atan(bh_a / r);
      double sph = std::sin(ph);
      double cph = std::cos(ph);

      // Calculate metric
      double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
      double sigma = r * r + bh_a * bh_a * cth * cth;
      double gtt_bl = -(1.0 + 2.0 * bh_m * r * (r * r + bh_a * bh_a) / (delta * sigma));
      double gtph_bl = -2.0 * bh_m * bh_a * r / (delta * sigma);
      double grr_bl = delta / sigma;
      double gthth_bl = 1.0 / sigma;
      double gphph_bl = (sigma - 2.0 * bh_m * r) / (delta * sigma * sth * sth);

      // Calculate angular momentum (C 6)
      double ll = formula_l0 / (1.0 + rr) * std::pow(rr, 1.0 + formula_q);

      // Calculate 4-velocity (C 7-8)
      double u_norm = 1.0 / std::sqrt(-gtt_bl + 2.0 * gtph_bl * ll - gphph_bl * ll * ll);
      double u_t_bl = -u_norm;
      double u_r_bl = 0.0;
      double u_th_bl = 0.0;
      double u_ph_bl = u_norm * ll;
      double ut_bl = gtt_bl * u_t_bl + gtph_bl * u_ph_bl;
      double ur_bl = grr_bl * u_r_bl;
      double uth_bl = gthth_bl * u_th_bl;
      double uph_bl = gtph_bl * u_t_bl + gphph_bl * u_ph_bl;
      double ut = ut_bl + 2.0 * bh_m * r / delta * ur_bl;
      double ur = ur_bl;
      double uth = uth_bl;
      double uph = uph_bl + bh_a / delta * ur_bl;
      double u0 = ut;
      double u1 =
          sth * cph * ur + cth * (r * cph - bh_a * sph) * uth + sth * (-r * sph - bh_a * cph) * uph;
      double u2 =
          sth * sph * ur + cth * (r * sph + bh_a * cph) * uth + sth * (r * cph - bh_a * sph) * uph;
      double u3 = cth * ur - r * sth * uth;

      // Calculate fluid-frame number density (C 5)
      double n_n0_fluid =
          std::exp(-0.5 * (r * r / (formula_r0 * formula_r0) + formula_h * formula_h * cth * cth));

      // Calculate frequency in CGS units
      double nu_fluid_cgs = -(u0 * k_0 + u1 * k_1 + u2 * k_2 + u3 * k_3) * momentum_factor;

      // Calculate emission coefficient in CGS units (C 9-10)
      double j_nu_fluid_cgs =
          formula_cn0 * n_n0_fluid * std::pow(nu_fluid_cgs / formula_nup, -formula_alpha);
      j_i_local(m,n) = j_nu_fluid_cgs / (nu_fluid_cgs * nu_fluid_cgs);

      // Calculate absorption coefficient in CGS units (C 11-12)
      double alpha_nu_fluid_cgs = formula_a * formula_cn0 * n_n0_fluid
          * std::pow(nu_fluid_cgs / formula_nup, -formula_beta - formula_alpha);
      alpha_i_local(m,n) = alpha_nu_fluid_cgs * nu_fluid_cgs;
    }
  }
  return;
}
