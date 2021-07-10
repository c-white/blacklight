// Blacklight ray tracer - simulation sampling and radiation integration

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // cbrt, cos, cosh, cyl_bessel_k, exp, expm1, pow, sin, sinh, sqrt, tanh
#include <complex>    // complex
#include <limits>     // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "ray_tracer.hpp"
#include "../blacklight.hpp"        // math, physics, enums
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function for resampling simulation cell data onto rays.
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes image_steps, sample_flags, sample_num, sample_pos, and sample_len have been set.
//   Allocates and initializes sample_rho, sample_pgas or sample_kappa, sample_uu1, sample_uu2,
//       sample_uu3, sample_bb1, sample_bb2, and sample_bb3.
//   If simulation_interp == false, uses primitives from cell containing geodesic sample point.
//   If simulation_interp == true and simulation_block_interp == false, performs trilinear
//       interpolation to geodesic sample point from cell centers, using only data within the same
//       block of cells (i.e. sometimes using extrapolation near block edges).
//   If simulation_interp == true and simulation_block_interp == true, performs trilinear
//       interpolation after obtaining anchor points possibly from neighboring blocks, even at
//       different refinement levels, or across the periodic boundary in spherical coordinates.
void RayTracer::SampleSimulationAlongGeodesics()
{
  // Allocate resampling arrays
  sample_rho.Allocate(image_resolution, image_resolution, image_steps);
  if (plasma_model == PlasmaModel::ti_te_beta)
    sample_pgas.Allocate(image_resolution, image_resolution, image_steps);
  if (plasma_model == PlasmaModel::code_kappa)
    sample_kappa.Allocate(image_resolution, image_resolution, image_steps);
  sample_uu1.Allocate(image_resolution, image_resolution, image_steps);
  sample_uu2.Allocate(image_resolution, image_resolution, image_steps);
  sample_uu3.Allocate(image_resolution, image_resolution, image_steps);
  sample_bb1.Allocate(image_resolution, image_resolution, image_steps);
  sample_bb2.Allocate(image_resolution, image_resolution, image_steps);
  sample_bb3.Allocate(image_resolution, image_resolution, image_steps);

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
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        // Extract number of steps along this geodesic
        int num_steps = sample_num(m,l);

        // Set NaN fallback values if geodesic poorly terminated
        if (fallback_nan and sample_flags(m,l))
        {
          for (int n = 0; n < num_steps; n++)
          {
            sample_rho(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            if (plasma_model == PlasmaModel::ti_te_beta)
              sample_pgas(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            if (plasma_model == PlasmaModel::code_kappa)
              sample_kappa(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_uu1(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_uu2(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_uu3(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_bb1(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_bb2(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_bb3(m,l,n) = std::numeric_limits<float>::quiet_NaN();
          }
          continue;
        }

        // Go along geodesic
        for (int n = 0; n < num_steps; n++)
        {
          // Extract coordinates
          double x1 = sample_pos(m,l,n,1);
          double x2 = sample_pos(m,l,n,2);
          double x3 = sample_pos(m,l,n,3);

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
              if (x1 >= x1_min_temp and x1 <= x1_max_temp and x2 >= x2_min_temp
                  and x2 <= x2_max_temp and x3 >= x3_min_temp and x3 <= x3_max_temp)
                break;
            }

            // Set fallback values if off grid
            if (b_new == n_b)
            {
              if (fallback_nan)
              {
                sample_rho(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                if (plasma_model == PlasmaModel::ti_te_beta)
                  sample_pgas(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                if (plasma_model == PlasmaModel::code_kappa)
                  sample_kappa(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                sample_uu1(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                sample_uu2(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                sample_uu3(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                sample_bb1(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                sample_bb2(m,l,n) = std::numeric_limits<float>::quiet_NaN();
                sample_bb3(m,l,n) = std::numeric_limits<float>::quiet_NaN();
              }
              else
              {
                sample_rho(m,l,n) = fallback_rho;
                if (plasma_model == PlasmaModel::ti_te_beta)
                  sample_pgas(m,l,n) = fallback_pgas;
                if (plasma_model == PlasmaModel::code_kappa)
                  sample_kappa(m,l,n) = fallback_kappa;
                sample_uu1(m,l,n) = fallback_uu1;
                sample_uu2(m,l,n) = fallback_uu2;
                sample_uu3(m,l,n) = fallback_uu3;
                sample_bb1(m,l,n) = fallback_bb1;
                sample_bb2(m,l,n) = fallback_bb2;
                sample_bb3(m,l,n) = fallback_bb3;
              }
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

          // Resample values without interpolation
          if (not simulation_interp)
          {
            sample_rho(m,l,n) = grid_rho(b,k,j,i);
            if (plasma_model == PlasmaModel::ti_te_beta)
              sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
            if (plasma_model == PlasmaModel::code_kappa)
              sample_kappa(m,l,n) = grid_kappa(b,k,j,i);
            sample_uu1(m,l,n) = grid_uu1(b,k,j,i);
            sample_uu2(m,l,n) = grid_uu2(b,k,j,i);
            sample_uu3(m,l,n) = grid_uu3(b,k,j,i);
            sample_bb1(m,l,n) = grid_bb1(b,k,j,i);
            sample_bb2(m,l,n) = grid_bb2(b,k,j,i);
            sample_bb3(m,l,n) = grid_bb3(b,k,j,i);
          }

          // Resample values with intrablock interpolation
          if (simulation_interp and not simulation_block_interp)
          {
            // Calculate interpolation/extrapolation indices and coefficients
            int i_m = i == 0 or (i != n_i - 1 and x1 >= static_cast<double>(x1v(b,i))) ? i : i - 1;
            int j_m = j == 0 or (j != n_j - 1 and x2 >= static_cast<double>(x2v(b,j))) ? j : j - 1;
            int k_m = k == 0 or (k != n_k - 1 and x3 >= static_cast<double>(x3v(b,k))) ? k : k - 1;
            int i_p = i_m + 1;
            int j_p = j_m + 1;
            int k_p = k_m + 1;
            double f_i = (x1 - static_cast<double>(x1v(b,i_m)))
                / (static_cast<double>(x1v(b,i_p)) - static_cast<double>(x1v(b,i_m)));
            double f_j = (x2 - static_cast<double>(x2v(b,j_m)))
                / (static_cast<double>(x2v(b,j_p)) - static_cast<double>(x2v(b,j_m)));
            double f_k = (x3 - static_cast<double>(x3v(b,k_m)))
                / (static_cast<double>(x3v(b,k_p)) - static_cast<double>(x3v(b,k_m)));

            // Interpolate density
            double val_mmm = static_cast<double>(grid_rho(b,k_m,j_m,i_m));
            double val_mmp = static_cast<double>(grid_rho(b,k_m,j_m,i_p));
            double val_mpm = static_cast<double>(grid_rho(b,k_m,j_p,i_m));
            double val_mpp = static_cast<double>(grid_rho(b,k_m,j_p,i_p));
            double val_pmm = static_cast<double>(grid_rho(b,k_p,j_m,i_m));
            double val_pmp = static_cast<double>(grid_rho(b,k_p,j_m,i_p));
            double val_ppm = static_cast<double>(grid_rho(b,k_p,j_p,i_m));
            double val_ppp = static_cast<double>(grid_rho(b,k_p,j_p,i_p));
            sample_rho(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate gas pressure
            if (plasma_model == PlasmaModel::ti_te_beta)
            {
              val_mmm = static_cast<double>(grid_pgas(b,k_m,j_m,i_m));
              val_mmp = static_cast<double>(grid_pgas(b,k_m,j_m,i_p));
              val_mpm = static_cast<double>(grid_pgas(b,k_m,j_p,i_m));
              val_mpp = static_cast<double>(grid_pgas(b,k_m,j_p,i_p));
              val_pmm = static_cast<double>(grid_pgas(b,k_p,j_m,i_m));
              val_pmp = static_cast<double>(grid_pgas(b,k_p,j_m,i_p));
              val_ppm = static_cast<double>(grid_pgas(b,k_p,j_p,i_m));
              val_ppp = static_cast<double>(grid_pgas(b,k_p,j_p,i_p));
              sample_pgas(m,l,n) =
                  static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                  + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                  + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                  + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                  + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);
            }

            // Interpolate electron entropy
            if (plasma_model == PlasmaModel::code_kappa)
            {
              val_mmm = static_cast<double>(grid_kappa(b,k_m,j_m,i_m));
              val_mmp = static_cast<double>(grid_kappa(b,k_m,j_m,i_p));
              val_mpm = static_cast<double>(grid_kappa(b,k_m,j_p,i_m));
              val_mpp = static_cast<double>(grid_kappa(b,k_m,j_p,i_p));
              val_pmm = static_cast<double>(grid_kappa(b,k_p,j_m,i_m));
              val_pmp = static_cast<double>(grid_kappa(b,k_p,j_m,i_p));
              val_ppm = static_cast<double>(grid_kappa(b,k_p,j_p,i_m));
              val_ppp = static_cast<double>(grid_kappa(b,k_p,j_p,i_p));
              sample_kappa(m,l,n) =
                  static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                  + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                  + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                  + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                  + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);
            }

            // Interpolate x1-velocity
            val_mmm = static_cast<double>(grid_uu1(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_uu1(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_uu1(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_uu1(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_uu1(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_uu1(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_uu1(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_uu1(b,k_p,j_p,i_p));
            sample_uu1(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x2-velocity
            val_mmm = static_cast<double>(grid_uu2(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_uu2(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_uu2(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_uu2(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_uu2(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_uu2(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_uu2(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_uu2(b,k_p,j_p,i_p));
            sample_uu2(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x3-velocity
            val_mmm = static_cast<double>(grid_uu3(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_uu3(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_uu3(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_uu3(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_uu3(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_uu3(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_uu3(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_uu3(b,k_p,j_p,i_p));
            sample_uu3(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x1-field
            val_mmm = static_cast<double>(grid_bb1(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_bb1(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_bb1(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_bb1(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_bb1(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_bb1(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_bb1(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_bb1(b,k_p,j_p,i_p));
            sample_bb1(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x2-field
            val_mmm = static_cast<double>(grid_bb2(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_bb2(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_bb2(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_bb2(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_bb2(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_bb2(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_bb2(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_bb2(b,k_p,j_p,i_p));
            sample_bb2(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Interpolate x3-field
            val_mmm = static_cast<double>(grid_bb3(b,k_m,j_m,i_m));
            val_mmp = static_cast<double>(grid_bb3(b,k_m,j_m,i_p));
            val_mpm = static_cast<double>(grid_bb3(b,k_m,j_p,i_m));
            val_mpp = static_cast<double>(grid_bb3(b,k_m,j_p,i_p));
            val_pmm = static_cast<double>(grid_bb3(b,k_p,j_m,i_m));
            val_pmp = static_cast<double>(grid_bb3(b,k_p,j_m,i_p));
            val_ppm = static_cast<double>(grid_bb3(b,k_p,j_p,i_m));
            val_ppp = static_cast<double>(grid_bb3(b,k_p,j_p,i_p));
            sample_bb3(m,l,n) = static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * val_mmm
                + (1.0 - f_k) * (1.0 - f_j) * f_i * val_mmp
                + (1.0 - f_k) * f_j * (1.0 - f_i) * val_mpm + (1.0 - f_k) * f_j * f_i * val_mpp
                + f_k * (1.0 - f_j) * (1.0 - f_i) * val_pmm + f_k * (1.0 - f_j) * f_i * val_pmp
                + f_k * f_j * (1.0 - f_i) * val_ppm + f_k * f_j * f_i * val_ppp);

            // Account for possible invalid values
            if (sample_rho(m,l,n) <= 0.0f)
              sample_rho(m,l,n) = grid_rho(b,k,j,i);
            if (plasma_model == PlasmaModel::ti_te_beta and sample_pgas(m,l,n) <= 0.0f)
              sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
            if (plasma_model == PlasmaModel::code_kappa and sample_kappa(m,l,n) <= 0.0f)
              sample_kappa(m,l,n) = grid_kappa(b,k,j,i);
          }

          // Resample values with interblock interpolation
          if (simulation_interp and simulation_block_interp)
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
            double vals[8][8];
            FindNearbyVals(b, k_m, j_m, i_m, vals[0]);
            FindNearbyVals(b, k_m, j_m, i_p, vals[1]);
            FindNearbyVals(b, k_m, j_p, i_m, vals[2]);
            FindNearbyVals(b, k_m, j_p, i_p, vals[3]);
            FindNearbyVals(b, k_p, j_m, i_m, vals[4]);
            FindNearbyVals(b, k_p, j_m, i_p, vals[5]);
            FindNearbyVals(b, k_p, j_p, i_m, vals[6]);
            FindNearbyVals(b, k_p, j_p, i_p, vals[7]);

            // Perform interpolation
            sample_rho(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][0]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][0]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][0]
                + (1.0 - f_k) * f_j * f_i * vals[3][0]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][0]
                + f_k * (1.0 - f_j) * f_i * vals[5][0]
                + f_k * f_j * (1.0 - f_i) * vals[6][0] + f_k * f_j * f_i * vals[7][0]);
            if (plasma_model == PlasmaModel::ti_te_beta)
              sample_pgas(m,l,n) =
                  static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][1]
                  + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][1]
                  + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][1]
                  + (1.0 - f_k) * f_j * f_i * vals[3][1]
                  + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][1]
                  + f_k * (1.0 - f_j) * f_i * vals[5][1]
                  + f_k * f_j * (1.0 - f_i) * vals[6][1] + f_k * f_j * f_i * vals[7][1]);
            if (plasma_model == PlasmaModel::code_kappa)
              sample_kappa(m,l,n) =
                  static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][1]
                  + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][1]
                  + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][1]
                  + (1.0 - f_k) * f_j * f_i * vals[3][1]
                  + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][1]
                  + f_k * (1.0 - f_j) * f_i * vals[5][1]
                  + f_k * f_j * (1.0 - f_i) * vals[6][1] + f_k * f_j * f_i * vals[7][1]);
            sample_uu1(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][2]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][2]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][2]
                + (1.0 - f_k) * f_j * f_i * vals[3][2]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][2]
                + f_k * (1.0 - f_j) * f_i * vals[5][2]
                + f_k * f_j * (1.0 - f_i) * vals[6][2] + f_k * f_j * f_i * vals[7][2]);
            sample_uu2(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][3]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][3]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][3]
                + (1.0 - f_k) * f_j * f_i * vals[3][3]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][3]
                + f_k * (1.0 - f_j) * f_i * vals[5][3]
                + f_k * f_j * (1.0 - f_i) * vals[6][3] + f_k * f_j * f_i * vals[7][3]);
            sample_uu3(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][4]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][4]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][4]
                + (1.0 - f_k) * f_j * f_i * vals[3][4]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][4]
                + f_k * (1.0 - f_j) * f_i * vals[5][4]
                + f_k * f_j * (1.0 - f_i) * vals[6][4] + f_k * f_j * f_i * vals[7][4]);
            sample_bb1(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][5]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][5]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][5]
                + (1.0 - f_k) * f_j * f_i * vals[3][5]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][5]
                + f_k * (1.0 - f_j) * f_i * vals[5][5]
                + f_k * f_j * (1.0 - f_i) * vals[6][5] + f_k * f_j * f_i * vals[7][5]);
            sample_bb2(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][6]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][6]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][6]
                + (1.0 - f_k) * f_j * f_i * vals[3][6]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][6]
                + f_k * (1.0 - f_j) * f_i * vals[5][6]
                + f_k * f_j * (1.0 - f_i) * vals[6][6] + f_k * f_j * f_i * vals[7][6]);
            sample_bb3(m,l,n) =
                static_cast<float>((1.0 - f_k) * (1.0 - f_j) * (1.0 - f_i) * vals[0][7]
                + (1.0 - f_k) * (1.0 - f_j) * f_i * vals[1][7]
                + (1.0 - f_k) * f_j * (1.0 - f_i) * vals[2][7]
                + (1.0 - f_k) * f_j * f_i * vals[3][7]
                + f_k * (1.0 - f_j) * (1.0 - f_i) * vals[4][7]
                + f_k * (1.0 - f_j) * f_i * vals[5][7]
                + f_k * f_j * (1.0 - f_i) * vals[6][7] + f_k * f_j * f_i * vals[7][7]);

            // Account for possible invalid values
            if (sample_rho(m,l,n) <= 0.0f)
              sample_rho(m,l,n) = grid_rho(b,k,j,i);
            if (plasma_model == PlasmaModel::ti_te_beta and sample_pgas(m,l,n) <= 0.0f)
              sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
            if (plasma_model == PlasmaModel::code_kappa and sample_kappa(m,l,n) <= 0.0f)
              sample_kappa(m,l,n) = grid_kappa(b,k,j,i);
          }
        }
      }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for integrating unpolarized radiative transfer equation based on simulation data
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_dir, sample_len, sample_rho, sample_pgas or sample_kappa,
//       sample_uu1, sample_uu2, sample_uu3, sample_bb1, sample_bb2, and sample_bb3 have been set.
//   Allocates and initializes image.
//   Assumes x^0 is ignorable.
//   References beta-dependent temperature ratio electron model from 2016 AA 586 A38 (E1).
//   References entropy-based electron model from 2017 MNRAS 466 705 (E2).
//   References symphony paper 2016 ApJ 822 34 (S).
void RayTracer::IntegrateSimulationUnpolarizedRadiation()
{
  // Allocate image array
  image.Allocate(image_resolution, image_resolution);
  image.Zero();

  // Calculate units
  double x_unit = physics::gg_msun * simulation_m_msun / (physics::c * physics::c);
  double e_unit = simulation_rho_cgs * physics::c * physics::c;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        // Go through samples
        int num_steps = sample_num(m,l);
        for (int n = 0; n < num_steps; n++)
        {
          // Extract geodesic position and momentum
          double x1 = sample_pos(m,l,n,1);
          double x2 = sample_pos(m,l,n,2);
          double x3 = sample_pos(m,l,n,3);
          double p_0 = sample_dir(m,l,n,0);
          double p_1 = sample_dir(m,l,n,1);
          double p_2 = sample_dir(m,l,n,2);
          double p_3 = sample_dir(m,l,n,3);

          // Extract model variables
          double rho = sample_rho(m,l,n);
          double pgas = 0.0;
          if (plasma_model == PlasmaModel::ti_te_beta)
            pgas = sample_pgas(m,l,n);
          double kappa = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa = sample_kappa(m,l,n);
          double uu1 = sample_uu1(m,l,n);
          double uu2 = sample_uu2(m,l,n);
          double uu3 = sample_uu3(m,l,n);
          double bb1 = sample_bb1(m,l,n);
          double bb2 = sample_bb2(m,l,n);
          double bb3 = sample_bb3(m,l,n);

          // Skip contribution if magnetic field vanishes
          if (bb1 == 0.0 and bb2 == 0.0 and bb3 == 0.0)
            continue;

          // Calculate metric
          CovariantCoordinateMetric(x1, x2, x3, gcov);
          ContravariantCoordinateMetric(x1, x2, x3, gcon);

          // Calculate 4-velocity
          double gamma = std::sqrt(1.0 + gcov(1,1) * uu1 * uu1 + 2.0 * gcov(1,2) * uu1 * uu2
              + 2.0 * gcov(1,3) * uu1 * uu3 + gcov(2,2) * uu2 * uu2 + 2.0 * gcov(2,3) * uu2 * uu3
              + gcov(3,3) * uu3 * uu3);
          double alpha = 1.0 / std::sqrt(-gcon(0,0));
          double beta1 = -gcon(0,1) / gcon(0,0);
          double beta2 = -gcon(0,2) / gcon(0,0);
          double beta3 = -gcon(0,3) / gcon(0,0);
          double u0 = gamma / alpha;
          double u1 = uu1 - beta1 * gamma / alpha;
          double u2 = uu2 - beta2 * gamma / alpha;
          double u3 = uu3 - beta3 * gamma / alpha;
          double u_0 = gcov(0,0) * u0 + gcov(0,1) * u1 + gcov(0,2) * u2 + gcov(0,3) * u3;
          double u_1 = gcov(1,0) * u0 + gcov(1,1) * u1 + gcov(1,2) * u2 + gcov(1,3) * u3;
          double u_2 = gcov(2,0) * u0 + gcov(2,1) * u1 + gcov(2,2) * u2 + gcov(2,3) * u3;
          double u_3 = gcov(3,0) * u0 + gcov(3,1) * u1 + gcov(3,2) * u2 + gcov(3,3) * u3;

          // Calculate magnetic field strength
          double b0 = u_1 * bb1 + u_2 * bb2 + u_3 * bb3;
          double b1 = (bb1 + b0 * u1) / u0;
          double b2 = (bb2 + b0 * u2) / u0;
          double b3 = (bb3 + b0 * u3) / u0;
          double b_0 = gcov(0,0) * b0 + gcov(0,1) * b1 + gcov(0,2) * b2 + gcov(0,3) * b3;
          double b_1 = gcov(1,0) * b0 + gcov(1,1) * b1 + gcov(1,2) * b2 + gcov(1,3) * b3;
          double b_2 = gcov(2,0) * b0 + gcov(2,1) * b1 + gcov(2,2) * b2 + gcov(2,3) * b3;
          double b_3 = gcov(3,0) * b0 + gcov(3,1) * b1 + gcov(3,2) * b2 + gcov(3,3) * b3;
          double b_sq = b_0 * b0 + b_1 * b1 + b_2 * b2 + b_3 * b3;

          // Skip contribution if considered to be vacuum
          if (plasma_sigma_max >= 0.0 and b_sq / rho > plasma_sigma_max)
            continue;

          // Calculate fluid-frame photon momentum
          double p_1f = p_1 - u_1 / u_0 * p_0;
          double p_2f = p_2 - u_2 / u_0 * p_0;
          double p_3f = p_3 - u_3 / u_0 * p_0;
          double g1f1f = gcon(1,1) + u1 * u1;
          double g1f2f = gcon(1,2) + u1 * u2;
          double g1f3f = gcon(1,3) + u1 * u3;
          double g2f2f = gcon(2,2) + u2 * u2;
          double g2f3f = gcon(2,3) + u2 * u3;
          double g3f3f = gcon(3,3) + u3 * u3;
          double p1f = g1f1f * p_1f + g1f2f * p_2f + g1f3f * p_3f;
          double p2f = g1f2f * p_1f + g2f2f * p_2f + g2f3f * p_3f;
          double p3f = g1f3f * p_1f + g2f3f * p_2f + g3f3f * p_3f;

          // Calculate fluid-frame magnetic field
          double bb1f = b1;
          double bb2f = b2;
          double bb3f = b3;
          double g_1f1f =
              gcov(1,1) - 2.0 * u_1 / u_0 * gcov(0,1) + u_1 * u_1 / (u_0 * u_0) * gcov(0,0);
          double g_1f2f = gcov(1,2) - u_1 / u_0 * gcov(0,2) - u_2 / u_0 * gcov(0,1)
              + u_1 * u_2 / (u_0 * u_0) * gcov(0,0);
          double g_1f3f = gcov(1,3) - u_1 / u_0 * gcov(0,3) - u_3 / u_0 * gcov(0,1)
              + u_1 * u_3 / (u_0 * u_0) * gcov(0,0);
          double g_2f2f =
              gcov(2,2) - 2.0 * u_2 / u_0 * gcov(0,2) + u_2 * u_2 / (u_0 * u_0) * gcov(0,0);
          double g_2f3f = gcov(2,3) - u_2 / u_0 * gcov(0,3) - u_3 / u_0 * gcov(0,2)
              + u_2 * u_3 / (u_0 * u_0) * gcov(0,0);
          double g_3f3f =
              gcov(3,3) - 2.0 * u_3 / u_0 * gcov(0,3) + u_3 * u_3 / (u_0 * u_0) * gcov(0,0);
          double bb_1f = g_1f1f * bb1f + g_1f2f * bb2f + g_1f3f * bb3f;
          double bb_2f = g_1f2f * bb1f + g_2f2f * bb2f + g_2f3f * bb3f;
          double bb_3f = g_1f3f * bb1f + g_2f3f * bb2f + g_3f3f * bb3f;

          // Calculate fluid-frame angle between photon direction and magnetic field
          double pf_dot_pf = p_1f * p1f + p_2f * p2f + p_3f * p3f;
          double bbf_dot_bbf = bb_1f * bb1f + bb_2f * bb2f + bb_3f * bb3f;
          double pf_dot_bbf = p_1f * bb1f + p_2f * bb2f + p_3f * bb3f;
          double cos2_theta = std::min(pf_dot_bbf * pf_dot_bbf / (pf_dot_pf * bbf_dot_bbf), 1.0);
          double sin_theta = std::sqrt(1.0 - cos2_theta);

          // Calculate fluid-frame quantities in CGS units
          double nu_fluid_cgs = -(u0 * p_0 + u1 * p_1 + u2 * p_2 + u3 * p_3) * momentum_factor;
          double n_cgs = rho * simulation_rho_cgs / (plasma_mu * physics::m_p);
          double n_e_cgs = n_cgs / (1.0 + 1.0 / plasma_ne_ni);
          double bb_cgs = std::sqrt(4.0 * math::pi * b_sq * e_unit);
          double nu_c_cgs = physics::e * bb_cgs / (2.0 * math::pi * physics::m_e * physics::c);

          // Calculate electron temperature for model with T_i/T_e a function of beta (E1 1)
          double kb_tt_e_cgs = 0.0;
          double theta_e = 0.0;
          if (plasma_model == PlasmaModel::ti_te_beta)
          {
            double beta_inv = b_sq / (2.0 * pgas);
            double tt_rat = (plasma_rat_high + plasma_rat_low * beta_inv * beta_inv)
                / (1.0 + beta_inv * beta_inv);
            double kb_tt_tot_cgs = plasma_mu * physics::m_p * physics::c * physics::c * pgas / rho;
            kb_tt_e_cgs = (plasma_ne_ni + 1.0) / (plasma_ne_ni + tt_rat) * kb_tt_tot_cgs;
            theta_e = kb_tt_e_cgs / (physics::m_e * physics::c * physics::c);
          }

          // Calculate electron temperature for given electron entropy (E2 13)
          if (plasma_model == PlasmaModel::code_kappa)
          {
            double mu_e = plasma_mu * (1.0 + 1.0 / plasma_ne_ni);
            double rho_e = rho * physics::m_e / (mu_e * physics::m_p);
            double rho_kappa_e_cbrt = std::cbrt(rho_e * kappa);
            theta_e =
                1.0 / 5.0 * (std::sqrt(1.0 + 25.0 * rho_kappa_e_cbrt * rho_kappa_e_cbrt) - 1.0);
            kb_tt_e_cgs = theta_e * physics::m_e * physics::c * physics::c;
          }

          // Calculate emission coefficient in CGS units (S 24), skipping contribution if necessary
          double nu_s_cgs = 2.0 / 9.0 * nu_c_cgs * theta_e * theta_e * sin_theta;
          if (nu_s_cgs == 0.0)
            continue;
          double k2 = std::cyl_bessel_k(2.0, 1.0 / theta_e);
          if (k2 == 0.0)
            continue;
          double xx = nu_fluid_cgs / nu_s_cgs;
          double xx_factor = std::sqrt(xx) + std::sqrt(std::cbrt(32.0 * math::sqrt2 * xx));
          double j_nu_fluid_cgs = math::sqrt2 * math::pi * physics::e * physics::e * n_e_cgs
              * nu_s_cgs / (3.0 * k2 * physics::c) * xx_factor * xx_factor
              * std::exp(-std::cbrt(xx));
          double j_nu_invariant_cgs = j_nu_fluid_cgs / (nu_fluid_cgs * nu_fluid_cgs);

          // Calculate absorption coefficient in CGS units (S 25)
          double b_nu_cgs = 2.0 * physics::h * nu_fluid_cgs * nu_fluid_cgs * nu_fluid_cgs
              / (physics::c * physics::c) / std::expm1(physics::h * nu_fluid_cgs / kb_tt_e_cgs);
          double k_nu_fluid_cgs = j_nu_fluid_cgs / b_nu_cgs;
          double k_nu_invariant_cgs = k_nu_fluid_cgs * nu_fluid_cgs;

          // Calculate change in invariant intensity
          double delta_lambda = sample_len(m,l,n);
          double delta_lambda_cgs = delta_lambda * x_unit / momentum_factor;
          if (k_nu_invariant_cgs > 0.0)
          {
            double delta_tau_nu = k_nu_invariant_cgs * delta_lambda_cgs;
            double ss_nu_invariant_cgs = j_nu_invariant_cgs / k_nu_invariant_cgs;
            if (delta_tau_nu <= delta_tau_max)
              image(m,l) = std::exp(-delta_tau_nu)
                  * (image(m,l) + ss_nu_invariant_cgs * std::expm1(delta_tau_nu));
            else
              image(m,l) = ss_nu_invariant_cgs;
          }
          else
            image(m,l) += j_nu_invariant_cgs * delta_lambda_cgs;
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

//--------------------------------------------------------------------------------------------------

// Function for integrating polarized radiative transfer equation based on simulation data
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_dir, sample_len, sample_rho, sample_pgas or sample_kappa,
//       sample_uu1, sample_uu2, sample_uu3, sample_bb1, sample_bb2, and sample_bb3 have been set.
//   Allocates and initializes image.
//   Assumes x^0 is ignorable.
//   References beta-dependent temperature ratio electron model from 2016 AA 586 A38 (E1).
//   References entropy-based electron model from 2017 MNRAS 466 705 (E2).
//   References grtrans paper 2016 MNRAS 462 115 (G)
//   References symphony paper 2016 ApJ 822 34 (S).
//     J_V in (S 31) has an overall sign error that is corrected here and in the symphony code.
//   References ipole paper 2018 MNRAS 475 43 (I).
//   References 1985 SoPh 97 239 (L)
//   Emissivities, absorptivities, and rotativities are calculated in their invariant forms,
//       disagreeing with their definitions in (S) and (G) but agreeing with their usages in (I).
//   Uses the definitions of nu_c and nu_s from (S); nu_{s,alt} is what (G) calls nu_c, and nu_c is
//       what (G) calls nu_B.
//   Tetrad is chosen such that j_U, alpha_U, rho_U = 0.
//   Integration proceeds via Strang splitting as in (I).
void RayTracer::IntegrateSimulationPolarizedRadiation()
{
  // Allocate image array
  image.Allocate(4, image_resolution, image_resolution);
  image.Zero();

  // Allocate and initialize coherency tensor array
  Array<std::complex<double>> nn_con(image_resolution, image_resolution, 4, 4);
  nn_con.Zero();

  // Calculate units
  double x_unit = physics::gg_msun * simulation_m_msun / (physics::c * physics::c);
  double e_unit = simulation_rho_cgs * physics::c * physics::c;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    double delta_lambda_old;
    double kcon_old[4];
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);
    Array<double> connection(4, 4, 4);
    Array<double> connection_old(4, 4, 4);
    Array<double> tetrad(4, 4);
    Array<std::complex<double>> nn_con_temp(4, 4);
    Array<std::complex<double>> nn_tet_cov(4, 4);
    Array<std::complex<double>> nn_tet_con(4, 4);
    Array<double> jacobian(4, 4);

    // Go through pixels, integrating radiation
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        // Check number of steps
        int num_steps = sample_num(m,l);
        if (num_steps <= 0)
          continue;

        // Zero registers
        delta_lambda_old = 0.0;
        nn_con_temp.Zero();

        // Go through samples
        for (int n = 0; n < num_steps; n++)
        {
          // Initialize flags indicating how coupling to matter should proceed
          bool couple = true;
          bool full_couple = true;
          bool optically_thin = true;

          // Extract affine step size
          double delta_lambda = sample_len(m,l,n);
          double delta_lambda_new = delta_lambda;
          if (n < num_steps - 1)
            delta_lambda_new = sample_len(m,l,n+1);
          double delta_lambda_cgs = delta_lambda * x_unit / momentum_factor;

          // Extract geodesic position and covariant momentum
          double x1 = sample_pos(m,l,n,1);
          double x2 = sample_pos(m,l,n,2);
          double x3 = sample_pos(m,l,n,3);
          double kcov[4];
          kcov[0] = sample_dir(m,l,n,0);
          kcov[1] = sample_dir(m,l,n,1);
          kcov[2] = sample_dir(m,l,n,2);
          kcov[3] = sample_dir(m,l,n,3);

          // Extract model variables
          double rho = sample_rho(m,l,n);
          double pgas = 0.0;
          if (plasma_model == PlasmaModel::ti_te_beta)
            pgas = sample_pgas(m,l,n);
          double kappa = 0.0;
          if (plasma_model == PlasmaModel::code_kappa)
            kappa = sample_kappa(m,l,n);
          double uu1 = sample_uu1(m,l,n);
          double uu2 = sample_uu2(m,l,n);
          double uu3 = sample_uu3(m,l,n);
          double bb1 = sample_bb1(m,l,n);
          double bb2 = sample_bb2(m,l,n);
          double bb3 = sample_bb3(m,l,n);

          // Calculate metric and connection
          CovariantCoordinateMetric(x1, x2, x3, gcov);
          ContravariantCoordinateMetric(x1, x2, x3, gcon);
          CoordinateConnection(x1, x2, x3, connection);
          if (n == 0)
            for (int mu = 0; mu < 4; mu++)
              for (int alpha = 0; alpha < 4; alpha++)
                for (int beta = 0; beta < 4; beta++)
                  connection_old(mu,alpha,beta) = connection(mu,alpha,beta);
          else
            for (int mu = 0; mu < 4; mu++)
              for (int alpha = 0; alpha < 4; alpha++)
                for (int beta = 0; beta < 4; beta++)
                  connection_old(mu,alpha,beta) =
                      0.5 * (connection_old(mu,alpha,beta) + connection(mu,alpha,beta));

          // Calculate geodesic contravariant momentum
          double kcon[4] = {};
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              kcon[mu] += gcon(mu,nu) * kcov[nu];
          if (n == 0)
            for (int mu = 0; mu < 4; mu++)
              kcon_old[mu] = kcon[mu];
          else
            for (int mu = 0; mu < 4; mu++)
              kcon_old[mu] = 0.5 * (kcon_old[mu] + kcon[mu]);

          // Parallel-transport N by first half step
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
            {
              std::complex<double> dnn_dlambda = 0.0;
              for (int alpha = 0; alpha < 4; alpha++)
                for (int beta = 0; beta < 4; beta++)
                  dnn_dlambda -= kcon_old[alpha]
                      * (connection_old(mu,alpha,beta) * nn_con(m,l,beta,nu)
                      + connection_old(nu,alpha,beta) * nn_con(m,l,mu,beta));
              nn_con_temp(mu,nu) += dnn_dlambda * (delta_lambda_old + delta_lambda) / 2.0;
            }
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              nn_con(m,l,mu,nu) = nn_con_temp(mu,nu);

          // Skip coupling if magnetic field vanishes
          if (bb1 == 0.0 and bb2 == 0.0 and bb3 == 0.0)
          {
            couple = false;
            full_couple = false;
          }

          // Calculate contravariant and covariant velocity
          double uu0 = std::sqrt(1.0 + gcov(1,1) * uu1 * uu1 + 2.0 * gcov(1,2) * uu1 * uu2
              + 2.0 * gcov(1,3) * uu1 * uu3 + gcov(2,2) * uu2 * uu2 + 2.0 * gcov(2,3) * uu2 * uu3
              + gcov(3,3) * uu3 * uu3);
          double lapse = 1.0 / std::sqrt(-gcon(0,0));
          double shift1 = -gcon(0,1) / gcon(0,0);
          double shift2 = -gcon(0,2) / gcon(0,0);
          double shift3 = -gcon(0,3) / gcon(0,0);
          double ucon[4];
          ucon[0] = uu0 / lapse;
          ucon[1] = uu1 - shift1 * uu0 / lapse;
          ucon[2] = uu2 - shift2 * uu0 / lapse;
          ucon[3] = uu3 - shift3 * uu0 / lapse;
          double ucov[4] = {};
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              ucov[mu] += gcov(mu,nu) * ucon[nu];

          // Calculate contravariant and covariant magnetic field
          double bcon[4];
          bcon[0] = ucov[1] * bb1 + ucov[2] * bb2 + ucov[3] * bb3;
          bcon[1] = (bb1 + bcon[0] * ucon[1]) / ucon[0];
          bcon[2] = (bb2 + bcon[0] * ucon[2]) / ucon[0];
          bcon[3] = (bb3 + bcon[0] * ucon[3]) / ucon[0];
          double bcov[4] = {};
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              bcov[mu] += gcov(mu,nu) * bcon[nu];
          double b_sq = 0.0;
          for (int mu = 0; mu < 4; mu++)
            b_sq += bcov[mu] * bcon[mu];

          // Skip coupling if considered to be vacuum
          if (plasma_sigma_max >= 0.0 and b_sq / rho > plasma_sigma_max)
          {
            couple = false;
            full_couple = false;
          }

          // Calculate orthonormal tetrad
          double upcon[4] = {};
          if (bb1 == 0.0 and bb2 == 0.0 and bb3 == 0.0)
          {
            if (simulation_coord == Coordinates::sph_ks)
              upcon[2] = 1.0;
            if (simulation_coord == Coordinates::cart_ks)
              upcon[3] = 1.0;
          }
          else
            for (int mu = 0; mu < 4; mu++)
              upcon[mu] = bcon[mu];
          Tetrad(ucon, ucov, kcon, kcov, upcon, gcov, gcon, tetrad);

          // Transform N into orthonormal frame
          nn_tet_cov.Zero();
          for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
              for (int mu = 0; mu < 4; mu++)
                for (int nu = 0; nu < 4; nu++)
                  for (int alpha = 0; alpha < 4; alpha++)
                    for (int beta = 0; beta < 4; beta++)
                      nn_tet_cov(a,b) += tetrad(a,mu) * tetrad(b,nu) * gcov(mu,alpha)
                          * gcov(nu,beta) * nn_con(m,l,alpha,beta);

          // Calculate orthonormal-frame Stokes quantities before coupling to fluid (I 14)
          double ss_start[4];
          ss_start[0] = 0.5 * (nn_tet_cov(1,1) + nn_tet_cov(2,2)).real();
          ss_start[1] = 0.5 * (nn_tet_cov(1,1) - nn_tet_cov(2,2)).real();
          ss_start[2] = 0.5 * (nn_tet_cov(1,2) + nn_tet_cov(2,1)).real();
          ss_start[3] = 0.5 * (nn_tet_cov(2,1) - nn_tet_cov(1,2)).imag();

          // Calculate electron temperature for model with T_i/T_e a function of beta (E1 1)
          double kb_tt_e_cgs = 0.0;
          double theta_e = 0.0;
          if (couple and plasma_model == PlasmaModel::ti_te_beta)
          {
            double beta_inv = b_sq / (2.0 * pgas);
            double tt_rat = (plasma_rat_high + plasma_rat_low * beta_inv * beta_inv)
                / (1.0 + beta_inv * beta_inv);
            double kb_tt_tot_cgs = plasma_mu * physics::m_p * physics::c * physics::c * pgas / rho;
            kb_tt_e_cgs = (plasma_ne_ni + 1.0) / (plasma_ne_ni + tt_rat) * kb_tt_tot_cgs;
            theta_e = kb_tt_e_cgs / (physics::m_e * physics::c * physics::c);
          }

          // Calculate electron temperature for given electron entropy (E2 13)
          if (couple and plasma_model == PlasmaModel::code_kappa)
          {
            double mu_e = plasma_mu * (1.0 + 1.0 / plasma_ne_ni);
            double rho_e = rho * physics::m_e / (mu_e * physics::m_p);
            double rho_kappa_e_cbrt = std::cbrt(rho_e * kappa);
            theta_e =
                1.0 / 5.0 * (std::sqrt(1.0 + 25.0 * rho_kappa_e_cbrt * rho_kappa_e_cbrt) - 1.0);
            kb_tt_e_cgs = theta_e * physics::m_e * physics::c * physics::c;
          }

          // Calculate orthonormal-frame angle between wavevector and magnetic field
          double sin2_theta_b, sin_theta_b, cos_theta_b;
          if (couple)
          {
            double k_tet_1 = 0.0;
            double k_tet_2 = 0.0;
            double k_tet_3 = 0.0;
            double b_tet_1 = 0.0;
            double b_tet_2 = 0.0;
            double b_tet_3 = 0.0;
            for (int mu = 0; mu < 4; mu++)
            {
              k_tet_1 += tetrad(1,mu) * kcov[mu];
              k_tet_2 += tetrad(2,mu) * kcov[mu];
              k_tet_3 += tetrad(3,mu) * kcov[mu];
              b_tet_1 += tetrad(1,mu) * bcov[mu];
              b_tet_2 += tetrad(2,mu) * bcov[mu];
              b_tet_3 += tetrad(3,mu) * bcov[mu];
            }
            double k_sq_tet = k_tet_1 * k_tet_1 + k_tet_2 * k_tet_2 + k_tet_3 * k_tet_3;
            double b_sq_tet = b_tet_1 * b_tet_1 + b_tet_2 * b_tet_2 + b_tet_3 * b_tet_3;
            double k_b_tet = k_tet_1 * b_tet_1 + k_tet_2 * b_tet_2 + k_tet_3 * b_tet_3;
            double cos2_theta_b = std::min(k_b_tet * k_b_tet / (k_sq_tet * b_sq_tet), 1.0);
            sin2_theta_b = 1.0 - cos2_theta_b;
            sin_theta_b = std::sqrt(sin2_theta_b);
            cos_theta_b = std::sqrt(cos2_theta_b) * (k_b_tet >= 0.0 ? 1.0 : -1.0);
          }

          // Calculate other orthonormal-frame quantities
          double nu_cgs = 0.0;
          double nu_2_cgs, nu_3_cgs, n_e_cgs, nu_c_cgs, nu_s_cgs, nu_s_alt_cgs;
          if (couple)
          {
            for (int mu = 0; mu < 4; mu++)
              nu_cgs -= kcov[mu] * ucon[mu];
            nu_cgs *= momentum_factor;
            nu_2_cgs = nu_cgs * nu_cgs;
            nu_3_cgs = nu_2_cgs * nu_cgs;
            double n_cgs = rho * simulation_rho_cgs / (plasma_mu * physics::m_p);
            n_e_cgs = n_cgs / (1.0 + 1.0 / plasma_ne_ni);
            double bb_cgs = std::sqrt(4.0 * math::pi * b_sq * e_unit);
            nu_c_cgs = physics::e * bb_cgs / (2.0 * math::pi * physics::m_e * physics::c);
            nu_s_cgs = 2.0 / 9.0 * nu_c_cgs * theta_e * theta_e * sin_theta_b;
            nu_s_alt_cgs = 27.0 / 4.0 * nu_s_cgs;
          }

          // Calculate thermal synchrotron emissivities (S 29,31)
          double j_i, j_q, j_u, j_v;
          if (couple)
          {
            double coefficient_j = n_e_cgs * physics::e * physics::e * nu_c_cgs / physics::c;
            double xx_j = nu_cgs / nu_s_cgs;
            double exp_neg_xx_1_3 = std::exp(-std::cbrt(xx_j));
            double xx_1_2 = std::sqrt(xx_j);
            double xx_1_6_times_2_11_12 = std::sqrt(std::cbrt(32.0 * math::sqrt2 * xx_j));
            double xx_9_25 = std::pow(xx_j, 9.0 / 25.0);
            double theta_e_24_25 = std::pow(theta_e, 24.0 / 25.0);
            double theta_e_3_5 = std::pow(theta_e, 3.0 / 5.0);
            double factor_i = xx_1_2 + xx_1_6_times_2_11_12;
            double factor_q = xx_1_2 + (7.0 * theta_e_24_25 + 35.0) / (10.0 * theta_e_24_25 + 75.0)
                * xx_1_6_times_2_11_12;
            double factor_v = 1.0 + (theta_e_3_5 / 25.0 + 7.0 / 10.0) * xx_9_25;
            double sin_theta_b_shift =
                sin_theta_b * std::cos(28.0 / 25.0) - cos_theta_b * std::sin(28.0 / 25.0);
            double coefficient_1 = math::sqrt2 * math::pi / 27.0;
            double coefficient_2 = (37.0 - 87.0 * sin_theta_b_shift) / (100.0 * (theta_e + 1.0));
            j_i = coefficient_j * exp_neg_xx_1_3 * coefficient_1 * sin_theta_b * factor_i * factor_i
                / nu_2_cgs;
            j_q = -coefficient_j * exp_neg_xx_1_3 * coefficient_1 * sin_theta_b * factor_q
                * factor_q / nu_2_cgs;
            j_u = 0.0;
            j_v = coefficient_j * exp_neg_xx_1_3 * coefficient_2 * std::pow(factor_v, 5.0 / 3.0)
                / nu_2_cgs;
          }

          // Calculate thermal synchrotron absorptivities from Kirchoff's law
          double alpha_i, alpha_q, alpha_u, alpha_v;
          if (couple)
          {
            double b_nu_cgs = 2.0 * physics::h * nu_cgs * nu_cgs * nu_cgs
                / (physics::c * physics::c) / std::expm1(physics::h * nu_cgs / kb_tt_e_cgs);
            alpha_i = nu_3_cgs * j_i / b_nu_cgs;
            alpha_q = nu_3_cgs * j_q / b_nu_cgs;
            alpha_u = 0.0;
            alpha_v = nu_3_cgs * j_v / b_nu_cgs;
            if (1.0 / (alpha_i * alpha_i) == std::numeric_limits<double>::infinity())
            {
              alpha_i = 0.0;
              alpha_q = 0.0;
              alpha_v = 0.0;
              full_couple = false;
            }
          }

          // Calculate optical depth
          double delta_tau;
          if (couple)
            delta_tau = alpha_i * delta_lambda_cgs;

          // Simplify coupling if optically thick
          if (couple)
            if (delta_tau > delta_tau_max)
              optically_thin = false;

          // Calculate thermal synchrotron rotativities
          double rho_q, rho_u, rho_v;
          if (couple)
          {
            double coefficient_rho = n_e_cgs * physics::e * physics::e * nu_c_cgs
                / (physics::m_e * physics::c);
            double kk_0 = std::cyl_bessel_k(0.0, 1.0 / theta_e);
            double kk_1 = std::cyl_bessel_k(1.0, 1.0 / theta_e);
            double kk_2 = std::cyl_bessel_k(2.0, 1.0 / theta_e);
            double xx_rho = 1.0
                / std::sqrt(3.0 / (2.0 * math::sqrt2) * 1.0e-3 * nu_cgs / nu_s_alt_cgs);  // (G B8)
            double xx_rho_exp = 0.011 * std::exp(-xx_rho / 47.2);
            double xx_rho_8_3 = 1.0e4 / std::cbrt(2.0 * std::sqrt(94143178827.0))
                * math::pi * std::pow(xx_rho, -8.0/3.0);
            double factor_f = 2.011 * std::exp(-std::pow(xx_rho, 1.035) / 4.7)
                - std::cos(xx_rho / 2.0) * std::exp(-std::pow(xx_rho, 1.2) / 2.73)
                - xx_rho_exp;                                                             // (G B6)
            double factor_fm = factor_f + (xx_rho_exp - xx_rho_8_3) * 0.5
                * (1.0 + std::tanh(10.0 * std::log(xx_rho / 120.0)));                     // (G B13)
            double delta_jj_5 =
                0.4379 * std::log(1.0 + 0.001858 * std::pow(xx_rho, 1.503));              // (G B15)
            rho_q = coefficient_rho * nu_c_cgs * sin2_theta_b / nu_2_cgs * factor_fm
                * (kk_1 / kk_2 + 6.0 * theta_e);                                          // (G B4)
            rho_u = 0.0;
            rho_v = 2.0 * coefficient_rho * cos_theta_b / nu_cgs * (kk_0 - delta_jj_5)
                / kk_2;                                                                   // (G B14)
          }

          // Calculate coefficients needed for coupling matrices
          double alpha_sq, alpha_p, rho_sq, rho_p, alpha_rho;
          double lambda_1, lambda_2, coefficient_theta, s;
          if (couple)
          {
            alpha_sq = alpha_q * alpha_q + alpha_u * alpha_u + alpha_v * alpha_v;
            alpha_p = std::sqrt(alpha_sq);
            rho_sq = rho_q * rho_q + rho_u * rho_u + rho_v * rho_v;
            rho_p = std::sqrt(rho_sq);
            alpha_rho = alpha_q * rho_q + alpha_u * rho_u + alpha_v * rho_v;
            double alpha_sq_rho_sq = alpha_sq - rho_sq;
            double lambda_a =
                std::sqrt(alpha_sq_rho_sq * alpha_sq_rho_sq / 4.0 + alpha_rho * alpha_rho);
            double lambda_b = alpha_sq_rho_sq / 2.0;
            lambda_1 = std::sqrt(lambda_a + lambda_b);
            lambda_2 = std::sqrt(lambda_a - lambda_b);
            coefficient_theta = lambda_1 * lambda_1 + lambda_2 * lambda_2;
            s = alpha_rho >= 0.0 ? 1.0 : -1.0;
          }

          // Calculate coupling matrix 1
          double mm_1[4][4] = {};
          if (full_couple)
            for (int a = 0; a < 4; a++)
              mm_1[a][a] = 1.0;

          // Calculate coupling matrix 2
          double mm_2[4][4] = {};
          if (full_couple)
          {
            mm_2[0][1] = lambda_2 * alpha_q - s * lambda_1 * rho_q;
            mm_2[0][2] = lambda_2 * alpha_u - s * lambda_1 * rho_u;
            mm_2[0][3] = lambda_2 * alpha_v - s * lambda_1 * rho_v;
            mm_2[1][2] = s * lambda_1 * alpha_v + lambda_2 * rho_v;
            mm_2[1][3] = -(s * lambda_1 * alpha_u + lambda_2 * rho_u);
            mm_2[1][2] = s * lambda_1 * alpha_q + lambda_2 * rho_q;
            mm_2[1][0] = mm_2[0][1];
            mm_2[2][0] = mm_2[0][2];
            mm_2[3][0] = mm_2[0][3];
            mm_2[2][1] = -mm_2[1][2];
            mm_2[3][1] = -mm_2[1][3];
            mm_2[3][2] = -mm_2[2][3];
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                mm_2[a][b] *= 1.0 / coefficient_theta;
          }

          // Calculate coupling matrix 3
          double mm_3[4][4] = {};
          if (full_couple)
          {
            mm_3[0][1] = lambda_1 * alpha_q + s * lambda_2 * rho_q;
            mm_3[0][2] = lambda_1 * alpha_u + s * lambda_2 * rho_u;
            mm_3[0][3] = lambda_1 * alpha_v + s * lambda_2 * rho_v;
            mm_3[1][2] = -(s * lambda_2 * alpha_v - lambda_1 * rho_v);
            mm_3[1][3] = s * lambda_2 * alpha_u - lambda_1 * rho_u;
            mm_3[1][2] = -(s * lambda_2 * alpha_q - lambda_1 * rho_q);
            mm_3[1][0] = mm_3[0][1];
            mm_3[2][0] = mm_3[0][2];
            mm_3[3][0] = mm_3[0][3];
            mm_3[2][1] = -mm_3[1][2];
            mm_3[3][1] = -mm_3[1][3];
            mm_3[3][2] = -mm_3[2][3];
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                mm_3[a][b] *= 1.0 / coefficient_theta;
          }

          // Calculate coupling matrix 4
          double mm_4[4][4] = {};
          if (full_couple)
          {
            mm_4[0][0] = (alpha_sq + rho_sq) / 2.0;
            mm_4[1][1] = alpha_q * alpha_q + rho_q * rho_q - (alpha_sq + rho_sq) / 2.0;
            mm_4[2][2] = alpha_u * alpha_u + rho_u * rho_u - (alpha_sq + rho_sq) / 2.0;
            mm_4[3][3] = alpha_v * alpha_v + rho_v * rho_v - (alpha_sq + rho_sq) / 2.0;
            mm_4[0][1] = alpha_v * rho_u - alpha_u * rho_v;
            mm_4[0][2] = alpha_q * rho_v - alpha_v * rho_q;
            mm_4[0][3] = alpha_u * rho_q - alpha_q * rho_u;
            mm_4[1][2] = alpha_q * alpha_u + rho_q * rho_u;
            mm_4[1][3] = alpha_v * alpha_q + rho_v * rho_q;
            mm_4[1][2] = alpha_u * alpha_v + rho_u * rho_v;
            mm_4[1][0] = -mm_4[0][1];
            mm_4[2][0] = -mm_4[0][2];
            mm_4[3][0] = -mm_4[0][3];
            mm_4[2][1] = mm_4[1][2];
            mm_4[3][1] = mm_4[1][3];
            mm_4[3][2] = mm_4[2][3];
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                mm_4[a][b] *= 2.0 / coefficient_theta;
          }

          // Calculate coupling polynomial O (L 10)
          double exp, sin, cos, sinh, cosh;
          double oo[4][4];
          if (full_couple and optically_thin)
          {
            exp = std::exp(-delta_tau);
            sin = std::sin(lambda_2 * delta_lambda_cgs);
            cos = std::cos(lambda_2 * delta_lambda_cgs);
            sinh = std::sinh(lambda_1 * delta_lambda_cgs);
            cosh = std::cosh(lambda_1 * delta_lambda_cgs);
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                oo[a][b] = exp * (0.5 * (mm_1[a][b] + mm_4[a][b]) * cosh
                    + 0.5 * (mm_1[a][b] - mm_4[a][b]) * cos - mm_2[a][b] * sin - mm_3[a][b] * sinh);
          }

          // Calculate coupling polynomial integral P (I 24)
          double pp[4][4];
          if (full_couple)
          {
            double f_1 = 1.0 / (alpha_i * alpha_i - lambda_1 * lambda_1);
            double f_2 = 1.0 / (alpha_i * alpha_i + lambda_2 * lambda_2);
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
              {
                double cosh_term =
                    -lambda_1 * f_1 * mm_3[a][b] + 0.5 * alpha_i * f_1 * (mm_1[a][b] + mm_4[a][b]);
                double cos_term =
                    -lambda_2 * f_2 * mm_2[a][b] + 0.5 * alpha_i * f_2 * (mm_1[a][b] - mm_4[a][b]);
                pp[a][b] = cosh_term + cos_term;
                if (optically_thin)
                {
                  double sin_term = -alpha_i * f_2 * mm_2[a][b]
                      - 0.5 * lambda_2 * f_2 * (mm_1[a][b] - mm_4[a][b]);
                  double sinh_term = -alpha_i * f_1 * mm_3[a][b]
                      + 0.5 * lambda_1 * f_1 * (mm_1[a][b] + mm_4[a][b]);
                  pp[a][b] -=
                      exp * (cosh_term * cosh + cos_term * cos + sin_term * sin + sinh_term * sinh);
                }
              }
          }

          // Couple matter to radiation over full step (I 23)
          double ss_end[4] = {};
          if (couple)
          {
            // Construct emissivity vector
            double j_s[4];
            j_s[0] = j_i;
            j_s[1] = j_q;
            j_s[2] = j_u;
            j_s[3] = j_v;

            // Construct absorptivity vector
            double alpha_s[4];
            alpha_s[0] = alpha_i;
            alpha_s[1] = alpha_q;
            alpha_s[2] = alpha_u;
            alpha_s[3] = alpha_v;

            // No absorptivity or rotativity
            if (alpha_i == 0.0 and rho_p == 0.0)
              for (int a = 0; a < 4; a++)
                ss_end[a] = ss_start[a] + j_s[a] * delta_lambda_cgs;

            // No polarized absorptivity or rotativity; nonzero absorptivity
            else if (alpha_p == 0.0 and rho_p == 0.0)
            {
              if (optically_thin)
              {
                double exp_neg = std::exp(-delta_tau);
                double expm1 = std::expm1(delta_tau);
                for (int a = 0; a < 4; a++)
                  ss_end[a] = exp_neg * (ss_start[a] + j_s[a] / alpha_i * expm1);
              }
              else
                for (int a = 0; a < 4; a++)
                  ss_end[a] = j_s[a] / alpha_i;
            }

            // No absorptivity; nonzero rotativity
            else if (alpha_i == 0.0)
            {
              double cos_rho = std::cos(rho_p * delta_lambda_cgs);
              double sin_rho = std::sin(rho_p * delta_lambda_cgs);
              double sin_sq_rho = std::sin(rho_p * delta_lambda_cgs / 2.0);
              sin_sq_rho = sin_sq_rho * sin_sq_rho;
              double rho_ss =
                  rho_q * ss_start[1] + rho_u * ss_start[2] + rho_v * ss_start[3];
              ss_end[0] = ss_start[0];                                              // (I A2)
              ss_end[1] = ss_start[1] * cos_rho
                  + 2.0 * rho_q * rho_ss / rho_sq * sin_sq_rho
                  + (rho_u * ss_start[3] - rho_v * ss_start[2]) * sin_rho;          // (I A3)
              ss_end[2] = ss_start[1] * cos_rho
                  + 2.0 * rho_u * rho_ss / rho_sq * sin_sq_rho
                  + (rho_v * ss_start[1] - rho_q * ss_start[3]) * sin_rho;          // (I A4)
              ss_end[3] = ss_start[3] * cos_rho
                  + 2.0 * rho_v * rho_ss / rho_sq * sin_sq_rho
                  + (rho_q * ss_start[2] - rho_u * ss_start[1]) * sin_rho;          // (I A5)
              for (int a = 0; a < 4; a++)
                ss_end[a] += j_s[a] * delta_lambda_cgs;
            }

            // No rotativity, nonzero polarized absorptivity
            else if (rho_p == 0.0)
            {
              if (optically_thin)
              {
                double exp_neg_i = std::exp(-delta_tau);
                double exp_neg_p = std::exp(-alpha_p * delta_lambda_cgs);
                double sinh_p = std::sinh(alpha_p * delta_lambda_cgs);
                double cosh_p = std::cosh(alpha_p * delta_lambda_cgs);
                double coshm1_p =
                    0.5 * (std::expm1(alpha_p * delta_lambda_cgs) + exp_neg_p - 1.0);
                double alpha_ss = alpha_q * ss_start[1] + alpha_u * ss_start[2]
                    + alpha_v * ss_start[3];
                double alpha_j = alpha_q * j_s[1] + alpha_u * j_s[2] + alpha_v * j_s[3];
                double alpha_i_p_factor = 1.0 / (alpha_i * alpha_i - alpha_sq);
                ss_end[0] = (ss_start[0] * cosh_p - alpha_ss / alpha_p * sinh_p)
                    * exp_neg_i + alpha_j * alpha_i_p_factor * (-1.0 + (alpha_i * sinh_p
                    + alpha_p * cosh_p) / alpha_p * exp_neg_p) + alpha_i * j_i
                    * alpha_i_p_factor * (1.0 - (alpha_i * cosh_p + alpha_p * sinh_p)
                    / alpha_i * exp_neg_p);                                               // (I A14)
                for (int a = 1; a < 4; a++)
                {
                  double term_1 = (ss_start[a] + alpha_s[a] * alpha_ss / alpha_sq
                      * coshm1_p - ss_start[0] * alpha_s[a] / alpha_p * sinh_p)
                      * exp_neg_i;
                  double term_2 = j_s[a] * (1.0 - exp_neg_i) / alpha_i;
                  double term_3 = alpha_j * alpha_s[a] / alpha_i * alpha_i_p_factor
                      * (1.0 - (1.0 - alpha_i * alpha_i / alpha_sq - alpha_i
                      / alpha_sq * (alpha_i * cosh_p + alpha_p * sinh_p))
                      * exp_neg_i);
                  double term_4 = j_i * alpha_s[a] / alpha_p * alpha_i_p_factor
                      * (-alpha_p + (alpha_p * cosh_p + alpha_i * sinh_p)
                      * exp_neg_i);
                  ss_end[a] = term_1 + term_2 + term_3 + term_4;                     // (I A15-A17)
                }
              }
              else
              {
                double alpha_j = alpha_q * j_s[1] + alpha_u * j_s[2] + alpha_v * j_s[3];
                ss_end[0] = (alpha_i * j_i - alpha_j) / (alpha_i * alpha_i - alpha_sq);
                for (int a = 1; a < 4; a++)
                  ss_end[a] = (j_s[a] - alpha_s[a] * ss_end[0]) / alpha_i;
              }
            }

            // Nonzero absorptivity and rotativity
            else
            {
              if (optically_thin)
                for (int a = 0; a < 4; a++)
                  for (int b = 0; b < 4; b++)
                    ss_end[a] += pp[a][b] * j_s[b] + oo[a][b] * ss_start[b];
              else
                for (int a = 0; a < 4; a++)
                  for (int b = 0; b < 4; b++)
                    ss_end[a] += pp[a][b] * j_s[b];
            }
          }

          // Account for no coupling
          else
            for (int a = 0; a < 4; a++)
              ss_end[a] = ss_start[a];

          // Calculate orthonormal-frame N after coupling to fluid (I 13)
          nn_tet_con.Zero();
          nn_tet_con(1,1) = ss_end[0] + ss_end[1];
          nn_tet_con(2,2) = ss_end[0] - ss_end[1];
          nn_tet_con(1,2) = ss_end[2] - math::i * ss_end[3];
          nn_tet_con(2,1) = ss_end[2] + math::i * ss_end[3];

          // Transform N into coordinate frame
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
            {
              nn_con(m,l,mu,nu) = 0.0;
              for (int a = 0; a < 4; a++)
                for (int b = 0; b < 4; b++)
                  nn_con(m,l,mu,nu) += tetrad(a,mu) * tetrad(b,nu) * nn_tet_con(a,b);
            }

          // Parallel-transport N by second half step
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
              nn_con_temp(mu,nu) = nn_con(m,l,mu,nu);
          for (int mu = 0; mu < 4; mu++)
            for (int nu = 0; nu < 4; nu++)
            {
              std::complex<double> dnn_dlambda = 0.0;
              for (int alpha = 0; alpha < 4; alpha++)
                for (int beta = 0; beta < 4; beta++)
                  dnn_dlambda -= kcon[alpha] * (connection(mu,alpha,beta) * nn_con_temp(beta,nu)
                      + connection(nu,alpha,beta) * nn_con_temp(mu,beta));
              nn_con(m,l,mu,nu) += dnn_dlambda * (delta_lambda + delta_lambda_new) / 4.0;
            }

          // Store values in registers for next step
          delta_lambda_old = delta_lambda;
          for (int mu = 0; mu < 4; mu++)
            kcon_old[mu] = kcon[mu];
          for (int mu = 0; mu < 4; mu++)
            for (int alpha = 0; alpha < 4; alpha++)
              for (int beta = 0; beta < 4; beta++)
                connection_old(mu,alpha,beta) = connection(mu,alpha,beta);
        }
      }

    // Go through pixels, transforming into camera frame
    #pragma omp for schedule(static)
    for (int m = 0; m < image_resolution; m++)
      for (int l = 0; l < image_resolution; l++)
      {
        // Extract geodesic position and covariant momentum
        double x = image_position(m,l,1);
        double y = image_position(m,l,2);
        double z = image_position(m,l,3);
        double kcov[4];
        kcov[0] = image_direction(m,l,0);
        kcov[1] = image_direction(m,l,1);
        kcov[2] = image_direction(m,l,2);
        kcov[3] = image_direction(m,l,3);

        // Calculate Jacobian of coordinate transform
        GeodesicCoordinateJacobian(x, y, z, jacobian);

        // Transform N into geodesic coordinates
        nn_con_temp.Zero();
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            for (int alpha = 0; alpha < 4; alpha++)
              for (int beta = 0; beta < 4; beta++)
                nn_con_temp(mu,nu) +=
                    jacobian(mu,alpha) * jacobian(nu,beta) * nn_con(m,l,alpha,beta);

        // Calculate metric and connection
        CovariantGeodesicMetric(x, y, z, gcov);
        ContravariantGeodesicMetric(x, y, z, gcon);

        // Calculate geodesic contravariant momentum
        double kcon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            kcon[mu] += gcon(mu,nu) * kcov[nu];

        // Calculate orientation
        double up_con[4];
        up_con[0] = camera_ucon[0] * camera_up_con_c[0] - (camera_ucov[1] * camera_up_con_c[1]
            + camera_ucov[2] * camera_up_con_c[2] + camera_ucov[3] * camera_up_con_c[3])
            / camera_ucov[0];
        up_con[1] = camera_up_con_c[1] + camera_ucon[1] * camera_up_con_c[0];
        up_con[2] = camera_up_con_c[2] + camera_ucon[2] * camera_up_con_c[0];
        up_con[3] = camera_up_con_c[3] + camera_ucon[3] * camera_up_con_c[0];

        // Calculate orthonormal tetrad
        Tetrad(camera_ucon, camera_ucov, kcon, kcov, up_con, gcov, gcon, tetrad);

        // Transform N into orthonormal frame
        nn_tet_cov.Zero();
        for (int a = 0; a < 4; a++)
          for (int b = 0; b < 4; b++)
            for (int mu = 0; mu < 4; mu++)
              for (int nu = 0; nu < 4; nu++)
                for (int alpha = 0; alpha < 4; alpha++)
                  for (int beta = 0; beta < 4; beta++)
                    nn_tet_cov(a,b) += tetrad(a,mu) * tetrad(b,nu) * gcov(mu,alpha) * gcov(nu,beta)
                        * nn_con_temp(alpha,beta);

        // Calculate orthonormal-frame Stokes quantities at camera location (I 14)
        image(0,m,l) = 0.5 * (nn_tet_cov(1,1) + nn_tet_cov(2,2)).real();
        image(1,m,l) = 0.5 * (nn_tet_cov(1,1) - nn_tet_cov(2,2)).real();
        image(2,m,l) = 0.5 * (nn_tet_cov(1,2) + nn_tet_cov(2,1)).real();
        image(3,m,l) = 0.5 * (nn_tet_cov(2,1) - nn_tet_cov(1,2)).imag();
      }

    // Transform invariant Stokes quantities (e.g. I_nu/nu^3) to standard ones (e.g. I_nu)
    double nu_cu = image_frequency * image_frequency * image_frequency;
    #pragma omp for schedule(static)
    for (int a = 0; a < 4; a++)
      for (int m = 0; m < image_resolution; m++)
        for (int l = 0; l < image_resolution; l++)
          image(a,m,l) *= nu_cu;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for finding cell values in a given or nearby block
// Inputs:
//   b: block index
//   k, j, i: cell indices for x^3, x^2, and x^1, possibly 1 beyond valid range
// Outputs:
//   vals: density, gas pressure or electron entropy, velocity, and magnetic field set
// Notes:
//   If requested cell is within block, values are copied from that cell.
//   If requested cell is in an adjacent block, values are copied from appropriate cell(s) there:
//     If adjacent block is at the same refinement level, the appropriate ghost cell is used.
//     If adjacent block is coarser, the coarse cell containing the appropriate ghost cell is used.
//     If adjacent block is finer, the equally-weighted average of the 8 cells comprising the
//         appropriate ghost cell is used.
//   If the requested cell is not on the grid, values are copied from the unique cell on the grid
//       closest to the appropriate ghost cell, effectively resulting in constant (rather than
//       linear) extrapolation near the edges of the grid.
//   In the case of simulation_coord == Coordinates::sph_ks, neighboring blocks are understood to
//       cross the periodic boundary in x^3 (phi), but the domain is not stitched together at the
//       poles.
void RayTracer::FindNearbyVals(int b, int k, int j, int i, double vals[8])
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
    vals[0] = grid_rho(b,k,j,i);
    if (plasma_model == PlasmaModel::ti_te_beta)
      vals[1] = grid_pgas(b,k,j,i);
    if (plasma_model == PlasmaModel::code_kappa)
      vals[1] = grid_kappa(b,k,j,i);
    vals[2] = grid_uu1(b,k,j,i);
    vals[3] = grid_uu2(b,k,j,i);
    vals[4] = grid_uu3(b,k,j,i);
    vals[5] = grid_bb1(b,k,j,i);
    vals[6] = grid_bb2(b,k,j,i);
    vals[7] = grid_bb3(b,k,j,i);
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
      vals[0] = grid_rho(b_alt,k_sought,j_sought,i_sought);
      if (plasma_model == PlasmaModel::ti_te_beta)
        vals[1] = grid_pgas(b_alt,k_sought,j_sought,i_sought);
      if (plasma_model == PlasmaModel::code_kappa)
        vals[1] = grid_kappa(b_alt,k_sought,j_sought,i_sought);
      vals[2] = grid_uu1(b_alt,k_sought,j_sought,i_sought);
      vals[3] = grid_uu2(b_alt,k_sought,j_sought,i_sought);
      vals[4] = grid_uu3(b_alt,k_sought,j_sought,i_sought);
      vals[5] = grid_bb1(b_alt,k_sought,j_sought,i_sought);
      vals[6] = grid_bb2(b_alt,k_sought,j_sought,i_sought);
      vals[7] = grid_bb3(b_alt,k_sought,j_sought,i_sought);
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
        vals[0] = grid_rho(b_alt,k_sought,j_sought,i_sought);
        if (plasma_model == PlasmaModel::ti_te_beta)
          vals[1] = grid_pgas(b_alt,k_sought,j_sought,i_sought);
        if (plasma_model == PlasmaModel::code_kappa)
          vals[1] = grid_kappa(b_alt,k_sought,j_sought,i_sought);
        vals[2] = grid_uu1(b_alt,k_sought,j_sought,i_sought);
        vals[3] = grid_uu2(b_alt,k_sought,j_sought,i_sought);
        vals[4] = grid_uu3(b_alt,k_sought,j_sought,i_sought);
        vals[5] = grid_bb1(b_alt,k_sought,j_sought,i_sought);
        vals[6] = grid_bb2(b_alt,k_sought,j_sought,i_sought);
        vals[7] = grid_bb3(b_alt,k_sought,j_sought,i_sought);
        return;
      }
  }

  // Find cells at finer level
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
      vals[0] = (static_cast<double>(grid_rho(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_rho(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_rho(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_rho(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_rho(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_rho(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_rho(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_rho(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      if (plasma_model == PlasmaModel::ti_te_beta)
        vals[1] = (static_cast<double>(grid_pgas(b_alt,k_sought,j_sought,i_sought))
            + static_cast<double>(grid_pgas(b_alt,k_sought,j_sought,i_sought+1))
            + static_cast<double>(grid_pgas(b_alt,k_sought,j_sought+1,i_sought))
            + static_cast<double>(grid_pgas(b_alt,k_sought,j_sought+1,i_sought+1))
            + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought,i_sought))
            + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought,i_sought+1))
            + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought+1,i_sought))
            + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      if (plasma_model == PlasmaModel::code_kappa)
        vals[1] = (static_cast<double>(grid_kappa(b_alt,k_sought,j_sought,i_sought))
            + static_cast<double>(grid_kappa(b_alt,k_sought,j_sought,i_sought+1))
            + static_cast<double>(grid_kappa(b_alt,k_sought,j_sought+1,i_sought))
            + static_cast<double>(grid_kappa(b_alt,k_sought,j_sought+1,i_sought+1))
            + static_cast<double>(grid_kappa(b_alt,k_sought+1,j_sought,i_sought))
            + static_cast<double>(grid_kappa(b_alt,k_sought+1,j_sought,i_sought+1))
            + static_cast<double>(grid_kappa(b_alt,k_sought+1,j_sought+1,i_sought))
            + static_cast<double>(grid_kappa(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      vals[2] = (static_cast<double>(grid_uu1(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_uu1(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_uu1(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_uu1(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_uu1(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_uu1(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_uu1(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_uu1(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      vals[3] = (static_cast<double>(grid_uu2(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_uu2(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_uu2(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_uu2(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_uu2(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_uu2(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_uu2(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_uu2(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      vals[4] = (static_cast<double>(grid_uu3(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_uu3(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_uu3(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_uu3(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_uu3(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_uu3(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_uu3(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_uu3(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      vals[5] = (static_cast<double>(grid_bb1(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_bb1(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_bb1(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_bb1(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_bb1(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_bb1(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_bb1(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_bb1(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      vals[6] = (static_cast<double>(grid_bb2(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_bb2(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_bb2(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_bb2(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_bb2(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_bb2(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_bb2(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_bb2(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      vals[7] = (static_cast<double>(grid_bb3(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_bb3(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_bb3(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_bb3(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_bb3(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_bb3(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_bb3(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_bb3(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
      return;
    }

  // Report grid inconsistency
  throw BlacklightException("Grid interpolation failed.");
  return;
}
