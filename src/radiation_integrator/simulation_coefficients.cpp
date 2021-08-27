// Blacklight radiation integrator - simulation radiative transfer coefficients

// C++ headers
#include <algorithm>  // min
#include <cmath>      // cbrt, cos, cyl_bessel_k, exp, expm1, pow, sin, sinh, sqrt, tanh
#include <limits>     // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // math, physics, enums
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for calculating radiative transfer coefficients along rays sampled from simulation
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Works on arrays appropriate for root level or adaptively refined regions.
//   Assumes geodesic_num_steps, sample_num (or sample_num_adaptive[adaptive_current_level]),
//       sample_pos (or sample_pos_adaptive[adaptive_current_level]), sample_dir (or
//       sample_dir_adaptive[adaptive_current_level]), sample_rho (or sample_rho_adaptive),
//       sample_pgas (or sample_pgas_adaptive) or sample_kappa (or sample_kappa_adaptive),
//       sample_uu1 (or sample_uu1_adaptive), sample_uu2 (or sample_uu2_adaptive), sample_uu3 (or
//       sample_uu3_adaptive), sample_bb1 (or sample_bb1_adaptive), sample_bb2 (or
//       sample_bb2_adaptive), and sample_bb3 (or sample_bb3_adaptive) have been set.
//   Allocates and initializes j_i (or j_i_adaptive) and alpha_i (or alpha_i_adaptive).
//   Allocates and initializes j_q (or j_q_adaptive), j_v (or j_v_adaptive), alpha_q (or
//       alpha_q_adaptive), alpha_v (or alpha_v_adaptive), rho_q (or rho_q_adaptive), and rho_v (or
//       rho_v_adaptive) if image_polarization == true.
//   References beta-dependent temperature ratio electron model from 2016 AA 586 A38 (E1).
//   References entropy-based electron model from 2017 MNRAS 466 705 (E2).
//   References symphony paper 2016 ApJ 822 34 (S).
//     J_V in (S 31) has an overall sign error that is corrected here and in the symphony code.
//   References grtrans paper 2016 MNRAS 462 115 (G)
//   Emissivities, absorptivities, and rotativities are calculated in their invariant forms,
//       disagreeing with their definitions in (S) and (G) but agreeing with their usages in 2018
//       MNRAS 475 43.
//   Uses the definitions of nu_c and nu_s from (S); nu_{s,alt} is what (G) calls nu_c, and nu_c is
//       what (G) calls nu_B.
//   Tetrad is chosen such that j_U, alpha_U, rho_U = 0.
//   Deallocates sample_rho_adaptive, sample_pgas_adaptive, and sample_kappa_adaptive.
//   Dealllocates sample_uu1_adaptive, sample_uu2_adaptive, sample_uu3_adaptive,
//       sample_bb1_adaptive, sample_bb2_adaptive, and sample_bb3_adaptive if
//       image_polarization == false.
void RadiationIntegrator::CalculateSimulationCoefficients()
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
    if (image_polarization)
    {
      j_q_adaptive.Allocate(num_pix, geodesic_num_steps_local);
      j_v_adaptive.Allocate(num_pix, geodesic_num_steps_local);
      alpha_q_adaptive.Allocate(num_pix, geodesic_num_steps_local);
      alpha_v_adaptive.Allocate(num_pix, geodesic_num_steps_local);
      rho_q_adaptive.Allocate(num_pix, geodesic_num_steps_local);
      rho_v_adaptive.Allocate(num_pix, geodesic_num_steps_local);
      j_q_adaptive.Zero();
      j_v_adaptive.Zero();
      alpha_q_adaptive.Zero();
      alpha_v_adaptive.Zero();
      rho_q_adaptive.Zero();
      rho_v_adaptive.Zero();
    }
  }
  else if (first_time)
  {
    j_i.Allocate(num_pix, geodesic_num_steps);
    alpha_i.Allocate(num_pix, geodesic_num_steps);
    j_i.Zero();
    alpha_i.Zero();
    if (image_polarization)
    {
      j_q.Allocate(num_pix, geodesic_num_steps);
      j_v.Allocate(num_pix, geodesic_num_steps);
      alpha_q.Allocate(num_pix, geodesic_num_steps);
      alpha_v.Allocate(num_pix, geodesic_num_steps);
      rho_q.Allocate(num_pix, geodesic_num_steps);
      rho_v.Allocate(num_pix, geodesic_num_steps);
      j_q.Zero();
      j_v.Zero();
      alpha_q.Zero();
      alpha_v.Zero();
      rho_q.Zero();
      rho_v.Zero();
    }
  }

  // Alias arrays
  Array<int> sample_num_local = sample_num;
  Array<double> sample_pos_local = sample_pos;
  Array<double> sample_dir_local = sample_dir;
  Array<float> sample_rho_local = sample_rho;
  Array<float> sample_pgas_local = sample_pgas;
  Array<float> sample_kappa_local = sample_kappa;
  Array<float> sample_uu1_local = sample_uu1;
  Array<float> sample_uu2_local = sample_uu2;
  Array<float> sample_uu3_local = sample_uu3;
  Array<float> sample_bb1_local = sample_bb1;
  Array<float> sample_bb2_local = sample_bb2;
  Array<float> sample_bb3_local = sample_bb3;
  Array<double> j_i_local = j_i;
  Array<double> j_q_local = j_q;
  Array<double> j_v_local = j_v;
  Array<double> alpha_i_local = alpha_i;
  Array<double> alpha_q_local = alpha_q;
  Array<double> alpha_v_local = alpha_v;
  Array<double> rho_q_local = rho_q;
  Array<double> rho_v_local = rho_v;
  if (adaptive_on and adaptive_current_level > 0)
  {
    sample_num_local = sample_num_adaptive[adaptive_current_level];
    sample_pos_local = sample_pos_adaptive[adaptive_current_level];
    sample_dir_local = sample_dir_adaptive[adaptive_current_level];
    sample_rho_local = sample_rho_adaptive;
    sample_pgas_local = sample_pgas_adaptive;
    sample_kappa_local = sample_kappa_adaptive;
    sample_uu1_local = sample_uu1_adaptive;
    sample_uu2_local = sample_uu2_adaptive;
    sample_uu3_local = sample_uu3_adaptive;
    sample_bb1_local = sample_bb1_adaptive;
    sample_bb2_local = sample_bb2_adaptive;
    sample_bb3_local = sample_bb3_adaptive;
    j_i_local = j_i_adaptive;
    j_q_local = j_q_adaptive;
    j_v_local = j_v_adaptive;
    alpha_i_local = alpha_i_adaptive;
    alpha_q_local = alpha_q_adaptive;
    alpha_v_local = alpha_v_adaptive;
    rho_q_local = rho_q_adaptive;
    rho_v_local = rho_v_adaptive;
  }

  // Calculate unit
  double e_unit = simulation_rho_cgs * physics::c * physics::c;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    Array<double> gcov_sim(4, 4);
    Array<double> gcon_sim(4, 4);
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);
    Array<double> jacobian(4, 4);
    Array<double> tetrad(4, 4);

    // Go through rays
    #pragma omp for schedule(static)
    for (int m = 0; m < num_pix; m++)
    {
      // Check number of steps
      int num_steps = sample_num_local(m);
      if (num_steps <= 0)
        continue;

      // Go through samples
      for (int n = 0; n < num_steps; n++)
      {
        // Extract geodesic position and covariant momentum
        double x1 = sample_pos_local(m,n,1);
        double x2 = sample_pos_local(m,n,2);
        double x3 = sample_pos_local(m,n,3);
        double kcov[4];
        kcov[0] = sample_dir_local(m,n,0);
        kcov[1] = sample_dir_local(m,n,1);
        kcov[2] = sample_dir_local(m,n,2);
        kcov[3] = sample_dir_local(m,n,3);

        // Extract model variables
        double rho = sample_rho_local(m,n);
        double pgas = 0.0;
        if (plasma_model == PlasmaModel::ti_te_beta)
          pgas = sample_pgas_local(m,n);
        double kappa = 0.0;
        if (plasma_model == PlasmaModel::code_kappa)
          kappa = sample_kappa_local(m,n);
        double uu1_sim = sample_uu1_local(m,n);
        double uu2_sim = sample_uu2_local(m,n);
        double uu3_sim = sample_uu3_local(m,n);
        double bb1_sim = sample_bb1_local(m,n);
        double bb2_sim = sample_bb2_local(m,n);
        double bb3_sim = sample_bb3_local(m,n);

        // Skip coupling if magnetic field vanishes
        if (bb1_sim == 0.0 and bb2_sim == 0.0 and bb3_sim == 0.0)
          continue;

        // Calculate simulation metric
        CovariantSimulationMetric(x1, x2, x3, gcov_sim);
        ContravariantSimulationMetric(x1, x2, x3, gcon_sim);

        // Calculate simulation velocity
        double uu0_sim = std::sqrt(1.0 + gcov_sim(1,1) * uu1_sim * uu1_sim
            + 2.0 * gcov_sim(1,2) * uu1_sim * uu2_sim + 2.0 * gcov_sim(1,3) * uu1_sim * uu3_sim
            + gcov_sim(2,2) * uu2_sim * uu2_sim + 2.0 * gcov_sim(2,3) * uu2_sim * uu3_sim
            + gcov_sim(3,3) * uu3_sim * uu3_sim);
        double lapse_sim = 1.0 / std::sqrt(-gcon_sim(0,0));
        double shift1_sim = -gcon_sim(0,1) / gcon_sim(0,0);
        double shift2_sim = -gcon_sim(0,2) / gcon_sim(0,0);
        double shift3_sim = -gcon_sim(0,3) / gcon_sim(0,0);
        double ucon_sim[4];
        ucon_sim[0] = uu0_sim / lapse_sim;
        ucon_sim[1] = uu1_sim - shift1_sim * uu0_sim / lapse_sim;
        ucon_sim[2] = uu2_sim - shift2_sim * uu0_sim / lapse_sim;
        ucon_sim[3] = uu3_sim - shift3_sim * uu0_sim / lapse_sim;
        double ucov_sim[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            ucov_sim[mu] += gcov_sim(mu,nu) * ucon_sim[nu];

        // Calculate simulation magnetic field
        double bcon_sim[4];
        bcon_sim[0] = ucov_sim[1] * bb1_sim + ucov_sim[2] * bb2_sim + ucov_sim[3] * bb3_sim;
        bcon_sim[1] = (bb1_sim + bcon_sim[0] * ucon_sim[1]) / ucon_sim[0];
        bcon_sim[2] = (bb2_sim + bcon_sim[0] * ucon_sim[2]) / ucon_sim[0];
        bcon_sim[3] = (bb3_sim + bcon_sim[0] * ucon_sim[3]) / ucon_sim[0];
        double bcov_sim[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            bcov_sim[mu] += gcov_sim(mu,nu) * bcon_sim[nu];
        double b_sq = 0.0;
        for (int mu = 0; mu < 4; mu++)
          b_sq += bcov_sim[mu] * bcon_sim[mu];

        // Skip coupling if considered to be vacuum
        if (plasma_sigma_max >= 0.0 and b_sq / rho > plasma_sigma_max)
          continue;

        // Calculate Jacobian of transformation from simulation to geodesic coordinates
        CoordinateJacobian(x1, x2, x3, jacobian);

        // Transform contravariant velocity and magnetic field to geodesic coordinates
        double ucon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            ucon[mu] += jacobian(mu,nu) * ucon_sim[nu];
        double bcon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            bcon[mu] += jacobian(mu,nu) * bcon_sim[nu];

        // Calculate geodesic metric
        CovariantGeodesicMetric(x1, x2, x3, gcov);
        ContravariantGeodesicMetric(x1, x2, x3, gcon);

        // Calculate geodesic contravariant momentum
        double kcon[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            kcon[mu] += gcon(mu,nu) * kcov[nu];

        // Calculate covariant velocity and magnetic field
        double ucov[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            ucov[mu] += gcov(mu,nu) * ucon[nu];
        double bcov[4] = {};
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            bcov[mu] += gcov(mu,nu) * bcon[nu];

        // Calculate orthonormal tetrad
        Tetrad(ucon, ucov, kcon, kcov, bcon, gcov, gcon, tetrad);

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
          theta_e = 1.0 / 5.0 * (std::sqrt(1.0 + 25.0 * rho_kappa_e_cbrt * rho_kappa_e_cbrt) - 1.0);
          kb_tt_e_cgs = theta_e * physics::m_e * physics::c * physics::c;
        }

        // Calculate orthonormal-frame angle between wavevector and magnetic field
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
        double sin2_theta_b = 1.0 - cos2_theta_b;
        double sin_theta_b = std::sqrt(sin2_theta_b);
        double cos_theta_b = std::sqrt(cos2_theta_b) * (k_b_tet >= 0.0 ? 1.0 : -1.0);

        // Calculate other orthonormal-frame quantities
        double nu_cgs = 0.0;
        for (int mu = 0; mu < 4; mu++)
          nu_cgs -= kcov[mu] * ucon[mu];
        nu_cgs *= momentum_factor;
        double nu_2_cgs = nu_cgs * nu_cgs;
        double nu_3_cgs = nu_2_cgs * nu_cgs;
        double n_cgs = rho * simulation_rho_cgs / (plasma_mu * physics::m_p);
        double n_e_cgs = n_cgs / (1.0 + 1.0 / plasma_ne_ni);
        double bb_cgs = std::sqrt(4.0 * math::pi * b_sq * e_unit);
        double nu_c_cgs = physics::e * bb_cgs / (2.0 * math::pi * physics::m_e * physics::c);
        double nu_s_cgs = 2.0 / 9.0 * nu_c_cgs * theta_e * theta_e * sin_theta_b;
        double nu_s_alt_cgs = 27.0 / 4.0 * nu_s_cgs;

        // Calculate thermal synchrotron emissivities (S 29,31)
        double coefficient_j = n_e_cgs * physics::e * physics::e * nu_c_cgs / physics::c;
        double xx_j = nu_cgs / nu_s_cgs;
        double exp_neg_xx_1_3 = std::exp(-std::cbrt(xx_j));
        double xx_1_2 = std::sqrt(xx_j);
        double xx_1_6_times_2_11_12 = std::sqrt(std::cbrt(32.0 * math::sqrt2 * xx_j));
        double factor_i = xx_1_2 + xx_1_6_times_2_11_12;
        double coefficient_1 = math::sqrt2 * math::pi / 27.0;
        j_i_local(m,n) = coefficient_j * exp_neg_xx_1_3 * coefficient_1 * sin_theta_b * factor_i
            * factor_i / nu_2_cgs;
        if (image_polarization)
        {
          double xx_9_25 = std::pow(xx_j, 9.0 / 25.0);
          double theta_e_24_25 = std::pow(theta_e, 24.0 / 25.0);
          double theta_e_3_5 = std::pow(theta_e, 3.0 / 5.0);
          double factor_q = xx_1_2 + (7.0 * theta_e_24_25 + 35.0) / (10.0 * theta_e_24_25 + 75.0)
              * xx_1_6_times_2_11_12;
          double factor_v = 1.0 + (theta_e_3_5 / 25.0 + 7.0 / 10.0) * xx_9_25;
          double sin_theta_b_shift =
              sin_theta_b * std::cos(28.0 / 25.0) - cos_theta_b * std::sin(28.0 / 25.0);
          double coefficient_2 = (37.0 - 87.0 * sin_theta_b_shift) / (100.0 * (theta_e + 1.0));
          j_q_local(m,n) = -coefficient_j * exp_neg_xx_1_3 * coefficient_1 * sin_theta_b * factor_q
              * factor_q / nu_2_cgs;
          j_v_local(m,n) = coefficient_j * exp_neg_xx_1_3 * coefficient_2
              * std::pow(factor_v, 5.0 / 3.0) / nu_2_cgs;
        }

        // Calculate thermal synchrotron absorptivities from Kirchoff's law
        double b_nu_cgs = 2.0 * physics::h * nu_cgs * nu_cgs * nu_cgs / (physics::c * physics::c)
            / std::expm1(physics::h * nu_cgs / kb_tt_e_cgs);
        alpha_i_local(m,n) = nu_3_cgs * j_i_local(m,n) / b_nu_cgs;
        if (image_polarization)
        {
          alpha_q_local(m,n) = nu_3_cgs * j_q_local(m,n) / b_nu_cgs;
          alpha_v_local(m,n) = nu_3_cgs * j_v_local(m,n) / b_nu_cgs;
        }

        // Account for numerical issues later arising from absorptivities being too small
        if (1.0 / (alpha_i_local(m,n) * alpha_i_local(m,n))
            == std::numeric_limits<double>::infinity())
        {
          alpha_i_local(m,n) = 0.0;
          if (image_polarization)
          {
            alpha_q_local(m,n) = 0.0;
            alpha_v_local(m,n) = 0.0;
          }
        }

        // Calculate thermal synchrotron rotativities (G B4,B6,B8,B13-B15)
        if (image_polarization)
        {
          double coefficient_rho =
              n_e_cgs * physics::e * physics::e * nu_c_cgs / (physics::m_e * physics::c);
          double kk_0 = std::cyl_bessel_k(0.0, 1.0 / theta_e);
          double kk_1 = std::cyl_bessel_k(1.0, 1.0 / theta_e);
          double kk_2 = std::cyl_bessel_k(2.0, 1.0 / theta_e);
          double xx_rho =
              1.0 / std::sqrt(3.0 / (2.0 * math::sqrt2) * 1.0e-3 * nu_cgs / nu_s_alt_cgs);
          double xx_rho_exp = 0.011 * std::exp(-xx_rho / 47.2);
          double xx_rho_8_3 = 1.0e4 / std::cbrt(2.0 * std::sqrt(94143178827.0)) * math::pi
              * std::pow(xx_rho, -8.0/3.0);
          double factor_f = 2.011 * std::exp(-std::pow(xx_rho, 1.035) / 4.7)
              - std::cos(xx_rho / 2.0) * std::exp(-std::pow(xx_rho, 1.2) / 2.73) - xx_rho_exp;
          double factor_fm = factor_f + (xx_rho_exp - xx_rho_8_3) * 0.5
              * (1.0 + std::tanh(10.0 * std::log(xx_rho / 120.0)));
          double delta_jj_5 = 0.4379 * std::log(1.0 + 0.001858 * std::pow(xx_rho, 1.503));
          rho_q_local(m,n) = coefficient_rho * nu_c_cgs * sin2_theta_b / nu_2_cgs * factor_fm
              * (kk_1 / kk_2 + 6.0 * theta_e);
          rho_v_local(m,n) =
              2.0 * coefficient_rho * cos_theta_b / nu_cgs * (kk_0 - delta_jj_5) / kk_2;
        }
      }
    }
  }

  // Free memory
  sample_rho_adaptive.Deallocate();
  sample_pgas_adaptive.Deallocate();
  sample_kappa_adaptive.Deallocate();
  if (not image_polarization)
  {
    sample_uu1_adaptive.Deallocate();
    sample_uu2_adaptive.Deallocate();
    sample_uu3_adaptive.Deallocate();
    sample_bb1_adaptive.Deallocate();
    sample_bb2_adaptive.Deallocate();
    sample_bb3_adaptive.Deallocate();
  }
  return;
}
