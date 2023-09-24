// Blacklight radiation integrator - simulation radiative transfer coefficients

// C++ headers
#include <algorithm>  // min
#include <cmath>      // cbrt, cos, cyl_bessel_k, exp, expm1, pow, sin, sinh, sqrt, tanh, tgamma
#include <limits>     // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Math, Physics, enums
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for calculating radiative transfer coefficients along rays sampled from simulation
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes geodesic_num_steps[adaptive_level], sample_num[adaptive_level],
//       sample_pos[adaptive_level], sample_dir[adaptive_level], sample_cut[adaptive_level],
//       sample_rho[adaptive_level], sample_pgas[adaptive_level], sample_kappa[adaptive_level] (if
//       needed), sample_uu1[adaptive_level], sample_uu2[adaptive_level],
//       sample_uu3[adaptive_level], sample_bb1[adaptive_level], sample_bb2[adaptive_level],
//       sample_bb3[adaptive_level], and momentum_factors[adaptive_level] have been set.
//   Allocates and initializes j_i[adaptive_level] if image_light == true or image_emission == true
//       or image_emission_ave == true.
//   Allocates and initializes alpha_i[adaptive_level] if image_light == true or image_tau == true
//       or image_tau_int == true.
//   Allocates and initializes j_q[adaptive_level], j_v[adaptive_level], alpha_q[adaptive_level],
//       alpha_v[adaptive_level], rho_q[adaptive_level], and rho_v[adaptive_level] if
//       image_light == true and image_polarization == true.
//   Allocates and initializes cell_values[adaptive_level] if image_lambda_ave == true
//       or image_emission_ave == true or image_tau_int == true or render_num_images > 0.
//   References beta-dependent temperature ratio electron model from 2016 AA 586 A38 (E1).
//   References entropy-based electron model from 2017 MNRAS 466 705 (E2).
//   References 2021 ApJ 921 17 (M) for transfer coefficients.
//   Transfer coefficients are calculated in their invariant forms, disagreeing with their
//       definitions in (M) but agreeing with their usages in 2018 MNRAS 475 43.
//   Faraday coefficients have additional trap for when Theta_e ~ 0, where rho_Q ~ 0 and rho_V is
//       given by cold-plasma rotation measure considerations, but numerically one might get NaN,
//       and the rho_V formula has the wrong asymptotic behavior.
//   Tetrad is chosen such that j_U, alpha_U, rho_U = 0.
//   Deallocates sample_cut[adaptive_level], sample_rho[adaptive_level],
//       sample_pgas[adaptive_level], and sample_kappa[adaptive_level] if adaptive_level > 0.
//   Dealllocates sample_uu1[adaptive_level], sample_uu2[adaptive_level],
//       sample_uu3[adaptive_level], sample_bb1[adaptive_level], sample_bb2[adaptive_level], and
//       sample_bb3[adaptive_level] if image_polarization == false and adaptive_level > 0.
void RadiationIntegrator::CalculateSimulationCoefficients()
{
  // Precalculate power-law values (M 38-42)
  if (first_time and plasma_power_frac != 0.0)
  {
    double var_a = std::pow(3.0, plasma_p / 2.0) * (plasma_p - 1.0);
    double var_b = 2.0 * (plasma_p + 1.0);
    double var_c =
        std::pow(plasma_gamma_min, 1.0 - plasma_p) - std::pow(plasma_gamma_max, 1.0 - plasma_p);
    double var_d = std::tgamma((3.0 * plasma_p - 1.0) / 12.0);
    double var_e = std::tgamma((3.0 * plasma_p + 19.0) / 12.0);
    double var_f = std::pow(3.0, (plasma_p + 1.0) / 2.0) * (plasma_p - 1.0) / 4.0;
    double var_g = std::tgamma((3.0 * plasma_p + 2.0) / 12.0);
    double var_h = std::tgamma((3.0 * plasma_p + 22.0) / 12.0);
    power_jj = var_a / var_b / var_c * var_d * var_e;
    power_aa = var_f / var_c * var_g * var_h;
    if (image_light and image_polarization)
    {
      double var_i = 2.0 * (plasma_p + 2.0) / (plasma_p + 1.0);
      double var_j = std::pow(plasma_gamma_min, -(plasma_p + 1.0));
      double var_k = std::log(plasma_gamma_min);
      power_jj_q = -(plasma_p + 1.0) / (plasma_p + 7.0 / 3.0);
      power_jj_v = 0.684 * std::pow(plasma_p, 0.49);
      power_aa_q = -std::pow(0.034 * plasma_p - 0.0344, 0.086);
      power_aa_v = std::pow(0.71 * plasma_p + 0.0352, 0.394);
      power_rho = (plasma_p - 1.0) / var_c;
      power_rho_q = -std::pow(plasma_gamma_min, 2.0 - plasma_p) / (plasma_p / 2.0 - 1.0);
      power_rho_v = var_i * var_j * var_k;
    }
  }

  // Precalculate kappa-distribution values (M 43,44,46-48,50-54)
  if (first_time and plasma_kappa_frac != 0.0)
  {
    double var_a = 4.0 * Math::pi * std::tgamma(plasma_kappa - 4.0 / 3.0);
    double var_b = std::pow(3.0, 7.0 / 3.0) * std::tgamma(plasma_kappa - 2.0);
    double var_c = std::pow(3.0, (plasma_kappa - 1.0) / 2.0);
    double var_d = (plasma_kappa - 2.0) * (plasma_kappa - 1.0) / 4.0;
    double var_e = std::tgamma(plasma_kappa / 4.0 - 1.0 / 3.0);
    double var_f = std::tgamma(plasma_kappa / 4.0 + 4.0 / 3.0);
    double var_g = std::pow(3.0, 1.0 / 6.0) * 10.0 / 41.0;
    double var_h = plasma_w * plasma_kappa;
    double var_i = 2.0 * Math::pi * std::pow(var_h, plasma_kappa - 10.0 / 3.0);
    double var_j = (plasma_kappa - 2.0) * (plasma_kappa - 1.0) * plasma_kappa;
    double var_k = 3.0 * plasma_kappa - 1.0;
    double var_l = std::tgamma(5.0 / 3.0);
    double var_m = Hypergeometric(plasma_kappa - 1.0 / 3.0, plasma_kappa + 1.0,
        plasma_kappa + 2.0 / 3.0, -var_h);
    double var_n = std::pow(Math::pi, 1.5) / 3.0;
    double var_o = var_j / (var_h * var_h * var_h);
    double var_p = 2.0 * std::tgamma(2.0 + plasma_kappa / 2.0) / (2.0 + plasma_kappa) - 1.0;
    kappa_jj_low = var_a / var_b;
    kappa_jj_high = var_c * var_d * var_e * var_f;
    kappa_jj_x_i = 3.0 * std::pow(plasma_kappa, -1.5);
    kappa_aa_low = var_g * var_i * var_j / var_k * var_l * var_m;
    kappa_aa_high = var_n * var_o * var_p;
    kappa_aa_x_i = std::pow(-1.75 + 1.6 * plasma_kappa, -0.86);
    if (image_light and image_polarization)
    {
      double var_q = 14.3 * std::pow(plasma_w, -0.928);
      double var_r = 169.0 * std::pow(plasma_kappa, -8.0) + 0.0052 * plasma_kappa - 0.0526
          + 47.0 / (200.0 * plasma_kappa);
      kappa_jj_low_q = 0.5;
      kappa_jj_low_v = 0.5625 * std::pow(plasma_kappa, -0.528) / plasma_w;
      kappa_jj_high_q = 0.64 + 0.02 * plasma_kappa;
      kappa_jj_high_v = 0.765625 * std::pow(plasma_kappa, -0.44) / plasma_w;
      kappa_jj_x_q = 3.7 * std::pow(plasma_kappa, -1.6);
      kappa_jj_x_v = kappa_jj_x_i;
      kappa_aa_low_q = 25.0 / 48.0;
      kappa_aa_low_v = 77.0 / (100.0 * plasma_w) * std::pow(plasma_kappa, -0.7);
      kappa_aa_high_i = std::pow(3.0 / plasma_kappa, 4.75) + 0.6;
      kappa_aa_high_q = 441.0 * std::pow(plasma_kappa, -5.76) + 0.55;
      kappa_aa_high_v = var_q * var_r;
      kappa_aa_x_q = 1.4 * std::pow(plasma_kappa, -1.15);
      kappa_aa_x_v = 1.22 * std::pow(plasma_kappa, -1.136) + 0.007;
      kappa_rho_v = std::cyl_bessel_k(0.0, 1.0 / plasma_w) / std::cyl_bessel_k(2.0, 1.0 / plasma_w);
      if (plasma_kappa < 4.0)
      {
        kappa_rho_frac = (plasma_kappa - 3.5) / (4.0 - 3.5);
        kappa_rho_q_low_a =
            17.0 * plasma_w + std::sqrt(plasma_w) * (-3.0 + 7.0 * std::exp(-5.0 * plasma_w));
        kappa_rho_q_low_b = -1.0 / 30.0;
        kappa_rho_q_low_c = 0.1;
        kappa_rho_q_low_d = -1.5;
        kappa_rho_q_low_e = 0.471;
        kappa_rho_q_high_a = 46.0 / 3.0 * plasma_w
            + std::sqrt(plasma_w) * (-5.0 / 3.0 + 17.0 / 3.0 * std::exp(-5.0 * plasma_w));
        kappa_rho_q_high_b = -1.0 / 18.0;
        kappa_rho_q_high_c = 1.0 / 6.0;
        kappa_rho_q_high_d = -1.75;
        kappa_rho_q_high_e = 0.5;
        kappa_rho_v_low_a = (plasma_w * plasma_w + 2.0 * plasma_w + 1.0)
            / (3.125 * plasma_w * plasma_w + 4.0 * plasma_w + 1.0);
        kappa_rho_v_low_b = 0.447;
        kappa_rho_v_high_a = (plasma_w * plasma_w + 54.0 * plasma_w + 50.0)
            / (30.0 / 11.0 * plasma_w * plasma_w + 134.0 * plasma_w + 50.0);
        kappa_rho_v_high_b = 0.391;
      }
      else if (plasma_kappa < 4.5)
      {
        kappa_rho_frac = (plasma_kappa - 4.0) / (4.5 - 4.0);
        kappa_rho_q_low_a = 46.0 / 3.0 * plasma_w
            + std::sqrt(plasma_w) * (-5.0 / 3.0 + 17.0 / 3.0 * std::exp(-5.0 * plasma_w));
        kappa_rho_q_low_b = -1.0 / 18.0;
        kappa_rho_q_low_c = 1.0 / 6.0;
        kappa_rho_q_low_d = -1.75;
        kappa_rho_q_low_e = 0.5;
        kappa_rho_q_high_a =
            14.0 * plasma_w + std::sqrt(plasma_w) * (-1.625 + 4.5 * std::exp(-5.0 * plasma_w));
        kappa_rho_q_high_b = -1.0 / 12.0;
        kappa_rho_q_high_c = 0.25;
        kappa_rho_q_high_d = -2.0;
        kappa_rho_q_high_e = 0.525;
        kappa_rho_v_low_a = (plasma_w * plasma_w + 54.0 * plasma_w + 50.0)
            / (30.0 / 11.0 * plasma_w * plasma_w + 134.0 * plasma_w + 50.0);
        kappa_rho_v_low_b = 0.391;
        kappa_rho_v_high_a = (plasma_w * plasma_w + 43.0 * plasma_w + 38.0)
            / (7.0 / 3.0 * plasma_w * plasma_w + 92.5 * plasma_w + 38.0);
        kappa_rho_v_high_b = 0.348;
      }
      else
      {
        kappa_rho_frac = (plasma_kappa - 4.5) / (5.0 - 4.5);
        kappa_rho_q_low_a =
            14.0 * plasma_w + std::sqrt(plasma_w) * (-1.625 + 4.5 * std::exp(-5.0 * plasma_w));
        kappa_rho_q_low_b = -1.0 / 12.0;
        kappa_rho_q_low_c = 0.25;
        kappa_rho_q_low_d = -2.0;
        kappa_rho_q_low_e = 0.525;
        kappa_rho_q_high_a =
            12.5 * plasma_w + std::sqrt(plasma_w) * (-1.0 + 5.0 * std::exp(-5.0 * plasma_w));
        kappa_rho_q_high_b = -0.125;
        kappa_rho_q_high_c = 0.375;
        kappa_rho_q_high_d = -2.25;
        kappa_rho_q_high_e = 0.541;
        kappa_rho_v_low_a = (plasma_w * plasma_w + 43.0 * plasma_w + 38.0)
            / (7.0 / 3.0 * plasma_w * plasma_w + 92.5 * plasma_w + 38.0);
        kappa_rho_v_low_b = 0.348;
        kappa_rho_v_high_a = (plasma_w + 13.0 / 14.0) / (2.0 * plasma_w + 13.0 / 14.0);
        kappa_rho_v_high_b = 0.313;
      }
    }
  }

  // Allocate arrays
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
  {
    if (image_light or image_emission or image_emission_ave)
      j_i[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
    if (image_light or image_tau or image_tau_int)
      alpha_i[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
    if (image_light and image_polarization)
    {
      j_q[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
      j_v[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
      alpha_q[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
      alpha_v[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
      rho_q[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
      rho_v[adaptive_level].Allocate(image_num_frequencies, num_pix,
          geodesic_num_steps[adaptive_level]);
    }
    if (image_lambda_ave or image_emission_ave or image_tau_int or render_num_images > 0)
      cell_values[adaptive_level].Allocate(CellValues::num_cell_values, num_pix,
          geodesic_num_steps[adaptive_level]);
  }
  j_i[adaptive_level].Zero();
  j_q[adaptive_level].Zero();
  j_v[adaptive_level].Zero();
  alpha_i[adaptive_level].Zero();
  alpha_q[adaptive_level].Zero();
  alpha_v[adaptive_level].Zero();
  rho_q[adaptive_level].Zero();
  rho_v[adaptive_level].Zero();
  cell_values[adaptive_level].SetNaN();

  // Calculate units
  double d_unit = simulation_rho_cgs;
  double e_unit = d_unit * Physics::c * Physics::c;
  double b_unit = std::sqrt(4.0 * Math::pi * e_unit);

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    double gcov_sim[4][4];
    double gcon_sim[4][4];
    double gcov[4][4];
    double gcon[4][4];
    double jacobian[4][4];
    double tetrad[4][4];

    // Go through rays and samples
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      int num_steps = sample_num[adaptive_level](m);
      for (int n = 0; n < num_steps; n++)
      {
        // Skip coupling if in cut region
        if (sample_cut[adaptive_level](m,n))
          continue;

        // Extract geodesic position and covariant momentum
        double x1 = sample_pos[adaptive_level](m,n,1);
        double x2 = sample_pos[adaptive_level](m,n,2);
        double x3 = sample_pos[adaptive_level](m,n,3);
        double kcov[4];
        kcov[0] = sample_dir[adaptive_level](m,n,0);
        kcov[1] = sample_dir[adaptive_level](m,n,1);
        kcov[2] = sample_dir[adaptive_level](m,n,2);
        kcov[3] = sample_dir[adaptive_level](m,n,3);

        // Extract model variables
        double rho = sample_rho[adaptive_level](m,n);
        double pgas = sample_pgas[adaptive_level](m,n);
        double kappa = 0.0;
        if (plasma_model == PlasmaModel::code_kappa)
          kappa = sample_kappa[adaptive_level](m,n);
        double uu1_sim = sample_uu1[adaptive_level](m,n);
        double uu2_sim = sample_uu2[adaptive_level](m,n);
        double uu3_sim = sample_uu3[adaptive_level](m,n);
        double bb1_sim = sample_bb1[adaptive_level](m,n);
        double bb2_sim = sample_bb2[adaptive_level](m,n);
        double bb3_sim = sample_bb3[adaptive_level](m,n);

        // Calculate densities and pressures
        double rho_cgs = rho * d_unit;
        double pgas_cgs = pgas * e_unit;
        double n_cgs = rho_cgs / (plasma_mu * Physics::m_p);
        double n_e_cgs = n_cgs / (1.0 + 1.0 / plasma_ne_ni);

        // Calculate simulation metric
        CovariantSimulationMetric(x1, x2, x3, gcov_sim);
        ContravariantSimulationMetric(x1, x2, x3, gcon_sim);

        // Calculate simulation velocity
        double uu0_sim = std::sqrt(1.0 + gcov_sim[1][1] * uu1_sim * uu1_sim
            + 2.0 * gcov_sim[1][2] * uu1_sim * uu2_sim + 2.0 * gcov_sim[1][3] * uu1_sim * uu3_sim
            + gcov_sim[2][2] * uu2_sim * uu2_sim + 2.0 * gcov_sim[2][3] * uu2_sim * uu3_sim
            + gcov_sim[3][3] * uu3_sim * uu3_sim);
        double lapse_sim = 1.0 / std::sqrt(-gcon_sim[0][0]);
        double shift1_sim = -gcon_sim[0][1] / gcon_sim[0][0];
        double shift2_sim = -gcon_sim[0][2] / gcon_sim[0][0];
        double shift3_sim = -gcon_sim[0][3] / gcon_sim[0][0];
        double ucon_sim[4];
        ucon_sim[0] = uu0_sim / lapse_sim;
        ucon_sim[1] = uu1_sim - shift1_sim * uu0_sim / lapse_sim;
        ucon_sim[2] = uu2_sim - shift2_sim * uu0_sim / lapse_sim;
        ucon_sim[3] = uu3_sim - shift3_sim * uu0_sim / lapse_sim;
        double ucov_sim[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            ucov_sim[mu] += gcov_sim[mu][nu] * ucon_sim[nu];

        // Calculate simulation magnetic field
        double bcon_sim[4];
        bcon_sim[0] = ucov_sim[1] * bb1_sim + ucov_sim[2] * bb2_sim + ucov_sim[3] * bb3_sim;
        bcon_sim[1] = (bb1_sim + bcon_sim[0] * ucon_sim[1]) / ucon_sim[0];
        bcon_sim[2] = (bb2_sim + bcon_sim[0] * ucon_sim[2]) / ucon_sim[0];
        bcon_sim[3] = (bb3_sim + bcon_sim[0] * ucon_sim[3]) / ucon_sim[0];
        double bcov_sim[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            bcov_sim[mu] += gcov_sim[mu][nu] * bcon_sim[nu];
        double b_sq = 0.0;
        for (int mu = 0; mu < 4; mu++)
          b_sq += bcov_sim[mu] * bcon_sim[mu];
        double bb_cgs = std::sqrt(b_sq) * b_unit;
        double sigma = b_sq / rho;
        double beta_inv = b_sq / (2.0 * pgas);

        // Calculate electron temperature for model with T_i/T_e a function of beta (E1 1)
        double kb_tt_e_cgs = std::numeric_limits<double>::quiet_NaN();
        double theta_e = std::numeric_limits<double>::quiet_NaN();
        if (plasma_thermal_frac != 0.0 and plasma_model == PlasmaModel::ti_te_beta)
        {
          double tt_rat = (plasma_rat_high + plasma_rat_low * beta_inv * beta_inv)
              / (1.0 + beta_inv * beta_inv);
          if (use_ipole_tpte)
          {
            double theta_e_unit = (adiabatic_gamma_elec - 1.0) * (adiabatic_gamma_ion - 1.0);
            theta_e_unit /= (adiabatic_gamma_ion - 1.0) + (adiabatic_gamma_elec - 1.0) * tt_rat;
            theta_e_unit *= Physics::m_p / Physics::m_e;
            theta_e = theta_e_unit * pgas / rho / (adiabatic_gamma - 1.0);
            kb_tt_e_cgs = theta_e * Physics::m_e * Physics::c * Physics::c;
          }
          else
          {
            double kb_tt_tot_cgs = plasma_mu * Physics::m_p * pgas_cgs / rho_cgs;
            kb_tt_e_cgs = (plasma_ne_ni + 1.0) / (plasma_ne_ni + tt_rat) * kb_tt_tot_cgs;
            theta_e = kb_tt_e_cgs / (Physics::m_e * Physics::c * Physics::c);
          }
        }

        // Calculate electron temperature for given electron entropy (E2 13)
        if (plasma_thermal_frac != 0.0 and plasma_model == PlasmaModel::code_kappa)
        {
          double mu_e = plasma_mu * (1.0 + 1.0 / plasma_ne_ni);
          double rho_e = rho * Physics::m_e / (mu_e * Physics::m_p);
          double rho_kappa_e_cbrt = std::cbrt(rho_e * kappa);
          theta_e = 1.0 / 5.0 * (std::sqrt(1.0 + 25.0 * rho_kappa_e_cbrt * rho_kappa_e_cbrt) - 1.0);
          kb_tt_e_cgs = theta_e * Physics::m_e * Physics::c * Physics::c;
        }

        // Skip coupling based on cell values
        if ((cut_rho_min >= 0.0 and rho_cgs < cut_rho_min)
            or (cut_rho_max >= 0.0 and rho_cgs > cut_rho_max)
            or (cut_n_e_min >= 0.0 and n_e_cgs < cut_n_e_min)
            or (cut_n_e_max >= 0.0 and n_e_cgs > cut_n_e_max)
            or (cut_p_gas_min >= 0.0 and pgas_cgs < cut_p_gas_min)
            or (cut_p_gas_max >= 0.0 and pgas_cgs > cut_p_gas_max)
            or (cut_theta_e_min >= 0.0 and theta_e < cut_theta_e_min)
            or (cut_theta_e_max >= 0.0 and theta_e > cut_theta_e_max)
            or (cut_b_min >= 0.0 and bb_cgs < cut_b_min)
            or (cut_b_max >= 0.0 and bb_cgs > cut_b_max)
            or (cut_sigma_min >= 0.0 and sigma < cut_sigma_min)
            or (cut_sigma_max >= 0.0 and sigma > cut_sigma_max)
            or (cut_beta_inverse_min >= 0.0 and beta_inv < cut_beta_inverse_min)
            or (cut_beta_inverse_max >= 0.0 and beta_inv > cut_beta_inverse_max))
          continue;

        // Record cell values
        if (image_lambda_ave or image_emission_ave or image_tau_int or render_num_images > 0)
        {
          cell_values[adaptive_level](static_cast<int>(CellValues::rho),m,n) = rho_cgs;
          cell_values[adaptive_level](static_cast<int>(CellValues::n_e),m,n) = n_e_cgs;
          cell_values[adaptive_level](static_cast<int>(CellValues::p_gas),m,n) = pgas_cgs;
          cell_values[adaptive_level](static_cast<int>(CellValues::theta_e),m,n) = theta_e;
          cell_values[adaptive_level](static_cast<int>(CellValues::bb),m,n) = bb_cgs;
          cell_values[adaptive_level](static_cast<int>(CellValues::sigma),m,n) = sigma;
          cell_values[adaptive_level](static_cast<int>(CellValues::beta_inv),m,n) = beta_inv;
        }

        // Skip remaining calculations if possible
        if (not (image_light or image_emission or image_tau or image_emission_ave or image_tau_int))
          continue;

        // Skip coupling if magnetic field vanishes
        if (bb1_sim == 0.0 and bb2_sim == 0.0 and bb3_sim == 0.0)
          continue;

        // Calculate Jacobian of transformation from simulation to geodesic coordinates
        CoordinateJacobian(x1, x2, x3, jacobian);

        // Transform contravariant velocity and magnetic field to geodesic coordinates
        double ucon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            ucon[mu] += jacobian[mu][nu] * ucon_sim[nu];
        double bcon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            bcon[mu] += jacobian[mu][nu] * bcon_sim[nu];

        // Calculate geodesic metric
        CovariantGeodesicMetric(x1, x2, x3, gcov);
        ContravariantGeodesicMetric(x1, x2, x3, gcon);

        // Calculate geodesic contravariant momentum
        double kcon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            kcon[mu] += gcon[mu][nu] * kcov[nu];

        // Calculate covariant velocity and magnetic field
        double ucov[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            ucov[mu] += gcov[mu][nu] * ucon[nu];
        double bcov[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            bcov[mu] += gcov[mu][nu] * bcon[nu];

        // Calculate orthonormal tetrad
        Tetrad(ucon, ucov, kcon, kcov, bcon, gcov, gcon, tetrad);

        // Calculate orthonormal-frame angle between wavevector and magnetic field
        double k_tet_1 = 0.0;
        double k_tet_2 = 0.0;
        double k_tet_3 = 0.0;
        double b_tet_1 = 0.0;
        double b_tet_2 = 0.0;
        double b_tet_3 = 0.0;
        for (int mu = 0; mu < 4; mu++)
        {
          k_tet_1 += tetrad[1][mu] * kcov[mu];
          k_tet_2 += tetrad[2][mu] * kcov[mu];
          k_tet_3 += tetrad[3][mu] * kcov[mu];
          b_tet_1 += tetrad[1][mu] * bcov[mu];
          b_tet_2 += tetrad[2][mu] * bcov[mu];
          b_tet_3 += tetrad[3][mu] * bcov[mu];
        }
        double k_sq_tet = k_tet_1 * k_tet_1 + k_tet_2 * k_tet_2 + k_tet_3 * k_tet_3;
        double b_sq_tet = b_tet_1 * b_tet_1 + b_tet_2 * b_tet_2 + b_tet_3 * b_tet_3;
        double k_b_tet = k_tet_1 * b_tet_1 + k_tet_2 * b_tet_2 + k_tet_3 * b_tet_3;
        double cos2_theta_b = std::min(k_b_tet * k_b_tet / (k_sq_tet * b_sq_tet), 1.0);
        double sin2_theta_b = 1.0 - cos2_theta_b;
        double sin_theta_b = std::sqrt(sin2_theta_b);
        double cos_theta_b = std::sqrt(cos2_theta_b) * (k_b_tet >= 0.0 ? 1.0 : -1.0);

        // Go through frequencies
        for (int l = 0; l < image_num_frequencies; l++)
        {
          // Calculate orthonormal-frame frequencies
          double nu_cgs = 0.0;
          for (int mu = 0; mu < 4; mu++)
            nu_cgs -= kcov[mu] * ucon[mu];
          nu_cgs *= image_frequencies(l) * momentum_factors[adaptive_level](m);
          double nu_2_cgs = nu_cgs * nu_cgs;
          double nu_c_cgs = Physics::e * bb_cgs / (2.0 * Math::pi * Physics::m_e * Physics::c);
          double nu_s_cgs = 2.0 / 9.0 * nu_c_cgs * theta_e * theta_e * sin_theta_b;

          // Calculate thermal synchrotron emissivities (M 28,30)
          double j_i_val;
          if (plasma_thermal_frac != 0.0)
          {
            double xx = nu_cgs / nu_s_cgs;
            double xx_1_2 = std::sqrt(xx);
            double xx_1_3 = std::cbrt(xx);
            double xx_1_6 = std::sqrt(xx_1_3);
            double coefficient = plasma_thermal_frac * n_e_cgs * Physics::e * Physics::e * nu_c_cgs
                / (Physics::c * nu_2_cgs) * std::exp(-xx_1_3);
            double var_a = Math::sqrt2 * Math::pi / 27.0 * sin_theta_b;
            double var_b = std::pow(2.0, 11.0 / 12.0);
            double var_c = xx_1_2 + var_b * xx_1_6;
            j_i_val = coefficient * var_a * var_c * var_c;
            if (image_light or image_emission or image_emission_ave)
              j_i[adaptive_level](l,m,n) = j_i_val;
            if (image_light and image_polarization)
            {
              double var_d = (7.0 * std::pow(theta_e, 0.96) + 35.0)
                  / (10.0 * std::pow(theta_e, 0.96) + 75.0) * var_b;
              double var_e = xx_1_2 + var_d * xx_1_6;
              double var_f = cos_theta_b / theta_e;
              double var_g = Math::pi / 3.0 + Math::pi / 3.0 * xx_1_3 + 2.0 / 300.0 * xx_1_2
                  + 2.0 / 19.0 * Math::pi * xx_1_3 * xx_1_3;
              j_q[adaptive_level](l,m,n) = -coefficient * var_a * var_e * var_e;
              j_v[adaptive_level](l,m,n) = coefficient * var_f * var_g;
            }
          }

          // Calculate thermal synchrotron absorptivities from Kirchoff's law (M 31)
          if (plasma_thermal_frac != 0.0)
          {
            // Calculate absorptivities
            double b_nu_nu_3_cgs = 2.0 * Physics::h / (Physics::c * Physics::c)
                / std::expm1(Physics::h * nu_cgs / kb_tt_e_cgs);
            if (image_light or image_tau or image_tau_int)
              alpha_i[adaptive_level](l,m,n) = j_i_val / b_nu_nu_3_cgs;
            if (image_light and image_polarization)
            {
              alpha_q[adaptive_level](l,m,n) = j_q[adaptive_level](l,m,n) / b_nu_nu_3_cgs;
              alpha_v[adaptive_level](l,m,n) = j_v[adaptive_level](l,m,n) / b_nu_nu_3_cgs;
            }

            // Account for numerical issues later arising from absorptivities being too small
            if ((image_light or image_tau or image_tau_int)
                and 1.0 / (alpha_i[adaptive_level](l,m,n) * alpha_i[adaptive_level](l,m,n))
                == std::numeric_limits<double>::infinity())
            {
              alpha_i[adaptive_level](l,m,n) = 0.0;
              if (image_light and image_polarization)
              {
                alpha_q[adaptive_level](l,m,n) = 0.0;
                alpha_v[adaptive_level](l,m,n) = 0.0;
              }
            }
          }

          // Calculate thermal synchrotron rotativities (M 33-37)
          if (plasma_thermal_frac != 0.0 and image_light and image_polarization)
          {
            double coefficient_q = -plasma_thermal_frac * n_e_cgs * Physics::e * Physics::e
                * nu_c_cgs * nu_c_cgs * sin2_theta_b / (Physics::m_e * Physics::c * nu_2_cgs);
            double coefficient_v = plasma_thermal_frac * 2.0 * n_e_cgs * Physics::e * Physics::e
                * nu_c_cgs * cos_theta_b / (Physics::m_e * Physics::c * nu_cgs);
            double factor_q = 0.0;
            double factor_v = 1.0;
            if (theta_e >= theta_e_zero)
            {
              double kk_0 = std::cyl_bessel_k(0.0, 1.0 / theta_e);
              double kk_1 = std::cyl_bessel_k(1.0, 1.0 / theta_e);
              double kk_2 = std::cyl_bessel_k(2.0, 1.0 / theta_e);
              double xx = nu_cgs / nu_s_cgs;
              double xx_neg_1_2 = 1.0 / std::sqrt(xx);
              double var_a = 2.011 * std::exp(-19.78 * std::pow(xx, -0.5175));
              double var_b = std::cos(39.89 * xx_neg_1_2) * std::exp(-70.16 * std::pow(xx, -0.6));
              double var_c = 0.011 * std::exp(-1.69 * xx_neg_1_2);
              double var_d = 0.003135 * std::pow(xx, 4.0 / 3.0);
              double var_e = 0.5 * (1.0 + std::tanh(10.0 * std::log(0.6648 * xx_neg_1_2)));
              double f_0 = var_a - var_b - var_c;
              double f_m = f_0 + (var_c - var_d) * var_e;
              double delta_jj_5 = 0.4379 * std::log(1.0 + 1.3414 * std::pow(xx, -0.7515));
              factor_q = f_m * (kk_1 / kk_2 + 6.0 * theta_e);
              factor_v = (kk_0 - delta_jj_5) / kk_2;
              factor_v = factor_v < 0.0 or factor_v > 1.0 ? 1.0 : factor_v;
            }
            rho_q[adaptive_level](l,m,n) = coefficient_q * factor_q;
            rho_v[adaptive_level](l,m,n) = coefficient_v * factor_v;
          }

          // Calculate power-law synchrotron emissivities (M 28,38)
          if (plasma_power_frac != 0.0 and (image_light or image_emission or image_emission_ave))
          {
            double var_a = std::pow(nu_cgs / (nu_c_cgs * sin_theta_b), -(plasma_p - 1.0) / 2.0);
            double coefficient = plasma_power_frac * n_e_cgs * Physics::e * Physics::e * nu_c_cgs
                / (Physics::c * nu_2_cgs) * power_jj * sin_theta_b * var_a;
            j_i[adaptive_level](l,m,n) += coefficient;
            if (image_light and image_polarization)
            {
              double var_b = cos_theta_b / sin_theta_b;
              double var_c = 1.0 / std::sqrt(nu_cgs / (3.0 * nu_c_cgs * sin_theta_b));
              j_q[adaptive_level](l,m,n) += coefficient * power_jj_q;
              j_v[adaptive_level](l,m,n) += coefficient * power_jj_v * var_b * var_c;
            }
          }

          // Calculate power-law synchrotron absorptivities (M 29,39)
          if (plasma_power_frac != 0.0 and (image_light or image_tau or image_tau_int))
          {
            double var_a = std::pow(nu_cgs / (nu_c_cgs * sin_theta_b), -(plasma_p + 2.0) / 2.0);
            double coefficient = plasma_power_frac * n_e_cgs * Physics::e * Physics::e
                / (Physics::m_e * Physics::c) * power_aa * var_a;
            alpha_i[adaptive_level](l,m,n) += coefficient;
            if (image_light and image_polarization)
            {
              double var_b = std::pow(3.1 * std::pow(sin_theta_b, -1.92) - 3.1, 0.512);
              double var_c = 1.0 / std::sqrt(nu_cgs / (nu_c_cgs * sin_theta_b));
              double var_d = cos_theta_b >= 0.0 ? 1.0 : -1.0;
              alpha_q[adaptive_level](l,m,n) += coefficient * power_aa_q;
              alpha_v[adaptive_level](l,m,n) += coefficient * power_aa_v * var_b * var_c * var_d;
            }
          }

          // Calculate power-law synchrotron rotativities (M 40-42)
          if (plasma_power_frac != 0.0 and image_light and image_polarization)
          {
            double var_a = n_e_cgs * Physics::e * Physics::e * nu_cgs
                / (Physics::m_e * Physics::c * nu_c_cgs * sin_theta_b);
            double var_b = nu_c_cgs * sin_theta_b / nu_cgs;
            double var_c = var_b * var_b;
            double var_d = var_c * var_b;
            double var_e = 1.0 - std::pow(2.0 * nu_c_cgs * plasma_gamma_min * plasma_gamma_min
                * sin_theta_b / (3.0 * nu_cgs), plasma_p / 2.0 - 1.0);
            double var_f = cos_theta_b / sin_theta_b;
            double coefficient = plasma_power_frac * power_rho * var_a;
            rho_q[adaptive_level](l,m,n) += coefficient * power_rho_q * var_d * var_e;
            rho_v[adaptive_level](l,m,n) += coefficient * power_rho_v * var_c * var_f;
          }

          // Calculate kappa-distribution synchrotron emissivities (M 28,43-46)
          if (plasma_kappa_frac != 0.0 and (image_light or image_emission or image_emission_ave))
          {
            double nu_kappa_cgs =
                nu_c_cgs * plasma_w * plasma_w * plasma_kappa * plasma_kappa * sin_theta_b;
            double xx = nu_cgs / nu_kappa_cgs;
            double var_a = plasma_kappa_frac * n_e_cgs * Physics::e * Physics::e * nu_c_cgs
                / (Physics::c * nu_2_cgs);
            double var_b = std::cbrt(xx) * sin_theta_b;
            double var_c = std::pow(xx, -(plasma_kappa - 2.0) / 2.0) * sin_theta_b;
            double coefficient_low = kappa_jj_low * var_a * var_b;
            double coefficient_high = kappa_jj_high * var_a * var_c;
            j_i[adaptive_level](l,m,n) += std::pow(std::pow(coefficient_low, -kappa_jj_x_i)
                + std::pow(coefficient_high, -kappa_jj_x_i), -1.0 / kappa_jj_x_i);
            if (image_light and image_polarization)
            {
              double var_d = std::pow(std::pow(sin_theta_b, -2.4) - 1.0, 0.48);
              double var_e = std::pow(xx, -0.35);
              double var_f = std::pow(std::pow(sin_theta_b, -2.5) - 1.0, 0.44);
              double var_g = 1.0 / std::sqrt(xx);
              double var_h = cos_theta_b >= 0.0 ? 1.0 : -1.0;
              double jj_q_low = coefficient_low * kappa_jj_low_q;
              double jj_v_low = coefficient_low * kappa_jj_low_v * var_d * var_e;
              double jj_q_high = coefficient_high * kappa_jj_high_q;
              double jj_v_high = coefficient_high * kappa_jj_high_v * var_f * var_g;
              j_q[adaptive_level](l,m,n) -= std::pow(std::pow(jj_q_low, -kappa_jj_x_q)
                  + std::pow(jj_q_high, -kappa_jj_x_q), -1.0 / kappa_jj_x_q);
              j_v[adaptive_level](l,m,n) += std::pow(std::pow(jj_v_low, -kappa_jj_x_v)
                  + std::pow(jj_v_high, -kappa_jj_x_v), -1.0 / kappa_jj_x_v) * var_h;
            }
          }

          // Calculate kappa-distribution synchrotron absoptivities (M 29,47-50)
          if (plasma_kappa_frac != 0.0 and (image_light or image_tau or image_tau_int))
          {
            double nu_kappa_cgs =
                nu_c_cgs * plasma_w * plasma_w * plasma_kappa * plasma_kappa * sin_theta_b;
            double xx = nu_cgs / nu_kappa_cgs;
            double var_a =
                plasma_kappa_frac * n_e_cgs * Physics::e * Physics::e / (Physics::m_e * Physics::c);
            double var_b = std::pow(xx, -2.0 / 3.0);
            double var_c = std::pow(xx, -(1.0 + plasma_kappa) / 2.0);
            double coefficient_low = kappa_aa_low * var_a * var_b;
            double coefficient_high = kappa_aa_high * var_a * var_c;
            double aa_i_low = coefficient_low;
            double aa_i_high = coefficient_high * kappa_aa_high_i;
            alpha_i[adaptive_level](l,m,n) += std::pow(std::pow(aa_i_low, -kappa_aa_x_i)
                + std::pow(aa_i_high, -kappa_aa_x_i), -1.0 / kappa_aa_x_i);
            if (image_light and image_polarization)
            {
              double var_d = std::pow(std::pow(sin_theta_b, -2.28) - 1.0, 0.446);
              double var_e = std::pow(xx, -0.35);
              double var_f = std::sqrt(std::pow(sin_theta_b, -2.05) - 1.0);
              double var_g = 1.0 / std::sqrt(xx);
              double var_h = cos_theta_b >= 0.0 ? 1.0 : -1.0;
              double aa_q_low = coefficient_low * kappa_aa_low_q;
              double aa_v_low = coefficient_low * kappa_aa_low_v * var_d * var_e;
              double aa_q_high = coefficient_high * kappa_aa_high_q;
              double aa_v_high = coefficient_high * kappa_aa_high_v * var_f * var_g;
              alpha_q[adaptive_level](l,m,n) -= std::pow(std::pow(aa_q_low, -kappa_aa_x_q)
                  + std::pow(aa_q_high, -kappa_aa_x_q), -1.0 / kappa_aa_x_q);
              alpha_v[adaptive_level](l,m,n) += std::pow(std::pow(aa_v_low, -kappa_aa_x_v)
                  + std::pow(aa_v_high, -kappa_aa_x_v), -1.0 / kappa_aa_x_v) * var_h;
            }
          }

          // Calculate kappa-distribution synchrotron rotativities (M 51-54)
          if (plasma_kappa_frac != 0.0 and image_light and image_polarization)
          {
            double nu_kappa_cgs =
                nu_c_cgs * plasma_w * plasma_w * plasma_kappa * plasma_kappa * sin_theta_b;
            double xx = nu_cgs / nu_kappa_cgs;
            double var_a = -plasma_kappa_frac * n_e_cgs * Physics::e * Physics::e * nu_c_cgs
                * nu_c_cgs * sin2_theta_b / (Physics::m_e * Physics::c * nu_2_cgs);
            double var_b = plasma_kappa_frac * 2.0 * n_e_cgs * Physics::e * Physics::e * nu_c_cgs
                * cos_theta_b / (Physics::m_e * Physics::c * nu_cgs);
            double var_c = 1.0 / std::sqrt(xx);
            double rho_q_low = var_a * kappa_rho_q_low_a * (1.0 - std::exp(kappa_rho_q_low_b
                * std::pow(xx, 0.84)) - std::sin(kappa_rho_q_low_c * xx)
                * std::exp(kappa_rho_q_low_d * std::pow(xx, kappa_rho_q_low_e)));
            double rho_q_high = var_a * kappa_rho_q_high_a * (1.0 - std::exp(kappa_rho_q_high_b
                * std::pow(xx, 0.84)) - std::sin(kappa_rho_q_high_c * xx)
                * std::exp(kappa_rho_q_high_d * std::pow(xx, kappa_rho_q_high_e)));
            double rho_v_low = kappa_rho_v * var_b * kappa_rho_v_low_a
                * (1.0 - 0.17 * std::log(1.0 + kappa_rho_v_low_b * var_c));
            double rho_v_high = kappa_rho_v * var_b * kappa_rho_v_high_a
                * (1.0 - 0.17 * std::log(1.0 + kappa_rho_v_high_b * var_c));
            rho_q[adaptive_level](l,m,n) +=
                (1.0 - kappa_rho_frac) * rho_q_low + kappa_rho_frac * rho_q_high;
            rho_v[adaptive_level](l,m,n) +=
                (1.0 - kappa_rho_frac) * rho_v_low + kappa_rho_frac * rho_v_high;
          }
        }
      }
    }
  }

  // Free memory
  if (adaptive_level > 0)
  {
    sample_cut[adaptive_level].Deallocate();
    sample_rho[adaptive_level].Deallocate();
    sample_pgas[adaptive_level].Deallocate();
    sample_kappa[adaptive_level].Deallocate();
    if (image_light and not image_polarization)
    {
      sample_uu1[adaptive_level].Deallocate();
      sample_uu2[adaptive_level].Deallocate();
      sample_uu3[adaptive_level].Deallocate();
      sample_bb1[adaptive_level].Deallocate();
      sample_bb2[adaptive_level].Deallocate();
      sample_bb3[adaptive_level].Deallocate();
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for evaluating hypergeometric function
// Inputs:
//   alpha, beta, gamma: function parameters
//   z: abscissa
// Outputs:
//   returned value: 2F1(alpha, beta, gamma, z)
// Notes:
//   Requires z < 0 and designed for use case of alpha = kappa - 1/3, beta = kappa + 1,
//       gamma = kappa + 2/3, and z = -w * kappa.
//   Evaluates Pfaff-transformation equivalent
//       (1 - z)^(-alpha) * 2F1(alpha, gamma - beta, gamma, z / (z - 1)).
//   Uses series expansion 2F1(a, b, c, x) = 1 + sum_{k=1}^infinity (a)_k (b)_k x^k / ((c)_k k!),
//       with (y)_k = y * (y + 1) * ... * (y + k - 1), that should converge to within 1% in the
//       range 3.5 <= kappa <= 5, 1 <= w <= 3.
double RadiationIntegrator::Hypergeometric(double alpha, double beta, double gamma, double z)
{
  // Parameters
  const int k_max = 10;

  // Perform transformation
  double a = alpha;
  double b = gamma - beta;
  double c = gamma;
  double x = z / (z - 1.0);

  // Prepare to accumulate result
  double result = 1.0;
  double a_k = 1.0;
  double b_k = 1.0;
  double c_k = 1.0;
  double xk = 1.0;
  double k_factorial = 1.0;

  // Sum series
  for (int k = 1; k <= k_max; k++)
  {
    a_k *= a + k - 1.0;
    b_k *= b + k - 1.0;
    c_k *= c + k - 1.0;
    xk *= x;
    k_factorial *= k;
    result += a_k * b_k * xk / (c_k * k_factorial);
  }

  // Scale result
  result *= std::pow(1.0 - z, -alpha);
  return result;
}
