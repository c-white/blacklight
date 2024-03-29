// Blacklight radiation integrator - formula radiative transfer coefficients

// C++ headers
#include <cmath>   // abs, acos, atan, atan2, cos, exp, pow, sin, sqrt
#include <limits>  // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Math
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on formula
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes sample_flags[adaptive_level], sample_num[adaptive_level], sample_pos[adaptive_level],
//       sample_dir[adaptive_level], and momentum_factors[adaptive_level] have been set.
//   Allocates and initializes j_i[adaptive_level] and alpha_i[adaptive_level].
//   References code comparison paper 2020 ApJ 897 148 (C).
void RadiationIntegrator::CalculateFormulaCoefficients()
{
  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
  {
    j_i[adaptive_level].Allocate(image_num_frequencies, num_pix,
        geodesic_num_steps[adaptive_level]);
    alpha_i[adaptive_level].Allocate(image_num_frequencies, num_pix,
        geodesic_num_steps[adaptive_level]);
  }
  j_i[adaptive_level].Zero();
  alpha_i[adaptive_level].Zero();

  // Go through rays in parallel
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < num_pix; m++)
  {
    // Check number of steps
    int num_steps = sample_num[adaptive_level](m);
    if (num_steps <= 0)
      continue;

    // Set pixel to NaN if ray has problem
    if (fallback_nan and sample_flags[adaptive_level](m))
    {
      for (int n = 0; n < num_steps; n++)
      {
        j_i[adaptive_level](m,n) = std::numeric_limits<double>::quiet_NaN();
        alpha_i[adaptive_level](m,n) = std::numeric_limits<double>::quiet_NaN();
      }
      continue;
    }

    // Go through samples
    for (int n = 0; n < num_steps; n++)
    {
      // Extract geodesic position and momentum
      double x = sample_pos[adaptive_level](m,n,1);
      double y = sample_pos[adaptive_level](m,n,2);
      double z = sample_pos[adaptive_level](m,n,3);
      double k_0 = sample_dir[adaptive_level](m,n,0);
      double k_1 = sample_dir[adaptive_level](m,n,1);
      double k_2 = sample_dir[adaptive_level](m,n,2);
      double k_3 = sample_dir[adaptive_level](m,n,3);

      // Cut outside camera radius
      double r = RadialGeodesicCoordinate(x, y, z);
      if (r > camera_r)
        continue;

      // Cut camera plane
      if (cut_omit_near or cut_omit_far)
      {
        double dot_product = x * camera_x[1] + y * camera_x[2] + z * camera_x[3];
        if ((cut_omit_near and dot_product > 0.0) or (cut_omit_far and dot_product < 0.0))
          continue;
      }

      // Cut spheres
      if ((cut_omit_in >= 0.0 and r < cut_omit_in) or (cut_omit_out >= 0.0 and r > cut_omit_out))
        continue;

      // Cut with respect to midplane
      if (cut_midplane_theta > 0.0 or cut_midplane_theta < 0.0)
      {
        double th = std::acos(z / r);
        if ((cut_midplane_theta > 0.0 and std::abs(th - Math::pi / 2.0) > cut_midplane_theta)
            or (cut_midplane_theta < 0.0 and std::abs(th - Math::pi / 2.0) < -cut_midplane_theta))
        {
          sample_cut[adaptive_level](m,n) = true;
          continue;
        }
      }
      if ((cut_midplane_z > 0.0 and std::abs(z) > cut_midplane_z)
          or (cut_midplane_z < 0.0 and std::abs(z) < -cut_midplane_z))
      {
        sample_cut[adaptive_level](m,n) = true;
        continue;
      }

      // Cut arbitrary plane
      if (cut_plane)
      {
        double dot_product = (x - cut_plane_origin_x) * cut_plane_normal_x
            + (y - cut_plane_origin_y) * cut_plane_normal_y
            + (z - cut_plane_origin_z) * cut_plane_normal_z;
        if (dot_product < 0.0)
          continue;
      }

      // Calculate curvilinear coordinates
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

      // Go through frequencies
      for (int l = 0; l < image_num_frequencies; l++)
      {
        // Calculate frequency in CGS units
        double nu_fluid_cgs = -(u0 * k_0 + u1 * k_1 + u2 * k_2 + u3 * k_3) * image_frequencies(l)
            * momentum_factors[adaptive_level](m);

        // Calculate emission coefficient in CGS units (C 9-10)
        double j_nu_fluid_cgs =
            formula_cn0 * n_n0_fluid * std::pow(nu_fluid_cgs / formula_nup, -formula_alpha);
        j_i[adaptive_level](l,m,n) = j_nu_fluid_cgs / (nu_fluid_cgs * nu_fluid_cgs);

        // Calculate absorption coefficient in CGS units (C 11-12)
        double alpha_nu_fluid_cgs = formula_a * formula_cn0 * n_n0_fluid
            * std::pow(nu_fluid_cgs / formula_nup, -formula_beta - formula_alpha);
        alpha_i[adaptive_level](l,m,n) = alpha_nu_fluid_cgs * nu_fluid_cgs;
      }
    }
  }
  return;
}
