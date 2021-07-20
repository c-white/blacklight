// Blacklight radiation integrator - polarized radiation integration

// C++ headers
#include <cmath>    // cos, cosh, exp, expm1, sin, sinh, sqrt
#include <complex>  // complex

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // math, physics
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for integrating polarized radiative transfer equation
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_len, j_i, j_q, j_v, alpha_i, alpha_q, alpha_v, rho_q, and rho_v have
//       been set.
//   Assumes sample_num, sample_pos, sample_dir, sample_len, sample_rho, sample_pgas or
//       sample_kappa, sample_uu1, sample_uu2, sample_uu3, sample_bb1, sample_bb2, and sample_bb3
//       have been set.
//   Allocates and initializes image.
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
void RadiationIntegrator::IntegratePolarizedRadiation()
{
  // Allocate image array
  image.Allocate(4, camera_num_pix);
  image.Zero();

  // Allocate and initialize coherency tensor array
  Array<std::complex<double>> nn_con(camera_num_pix, 4, 4);
  nn_con.Zero();

  // Calculate unit
  double x_unit = physics::gg_msun * mass_msun / (physics::c * physics::c);

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    double delta_lambda_old;
    double kcon_old[4];
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);
    Array<double> gcov_sim(4, 4);
    Array<double> gcon_sim(4, 4);
    Array<double> connection(4, 4, 4);
    Array<double> connection_old(4, 4, 4);
    Array<double> tetrad(4, 4);
    Array<std::complex<double>> nn_con_temp(4, 4);
    Array<std::complex<double>> nn_tet_cov(4, 4);
    Array<std::complex<double>> nn_tet_con(4, 4);
    Array<double> jacobian(4, 4);

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < camera_num_pix; m++)
    {
      // Check number of steps
      int num_steps = sample_num(m);
      if (num_steps <= 0)
        continue;

      // Zero registers
      delta_lambda_old = 0.0;
      nn_con_temp.Zero();

      // Go through samples
      for (int n = 0; n < num_steps; n++)
      {
        // Extract affine step size
        double delta_lambda = sample_len(m,n);
        double delta_lambda_new = delta_lambda;
        if (n < num_steps - 1)
          delta_lambda_new = sample_len(m,n+1);
        double delta_lambda_cgs = delta_lambda * x_unit / momentum_factor;

        // Extract geodesic position and covariant momentum
        double x1 = sample_pos(m,n,1);
        double x2 = sample_pos(m,n,2);
        double x3 = sample_pos(m,n,3);
        double kcov[4];
        kcov[0] = sample_dir(m,n,0);
        kcov[1] = sample_dir(m,n,1);
        kcov[2] = sample_dir(m,n,2);
        kcov[3] = sample_dir(m,n,3);

        // Extract model variables
        double uu1_sim = sample_uu1(m,n);
        double uu2_sim = sample_uu2(m,n);
        double uu3_sim = sample_uu3(m,n);
        double bb1_sim = sample_bb1(m,n);
        double bb2_sim = sample_bb2(m,n);
        double bb3_sim = sample_bb3(m,n);

        // Calculate geodesic metric and connection
        CovariantGeodesicMetric(x1, x2, x3, gcov);
        ContravariantGeodesicMetric(x1, x2, x3, gcon);
        GeodesicConnection(x1, x2, x3, connection);
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
                dnn_dlambda -= kcon_old[alpha] * (connection_old(mu,alpha,beta) * nn_con(m,beta,nu)
                    + connection_old(nu,alpha,beta) * nn_con(m,mu,beta));
            nn_con_temp(mu,nu) += dnn_dlambda * (delta_lambda_old + delta_lambda) / 2.0;
          }
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            nn_con(m,mu,nu) = nn_con_temp(mu,nu);

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
        double upcon[4] = {};
        if (bb1_sim == 0.0 and bb2_sim == 0.0 and bb3_sim == 0.0)
          upcon[3] = 1.0;
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
                    nn_tet_cov(a,b) += tetrad(a,mu) * tetrad(b,nu) * gcov(mu,alpha) * gcov(nu,beta)
                        * nn_con(m,alpha,beta);

        // Calculate orthonormal-frame Stokes quantities before coupling to fluid (I 14)
        double ss_start[4];
        ss_start[0] = 0.5 * (nn_tet_cov(1,1) + nn_tet_cov(2,2)).real();
        ss_start[1] = 0.5 * (nn_tet_cov(1,1) - nn_tet_cov(2,2)).real();
        ss_start[2] = 0.5 * (nn_tet_cov(1,2) + nn_tet_cov(2,1)).real();
        ss_start[3] = 0.5 * (nn_tet_cov(2,1) - nn_tet_cov(1,2)).imag();

        // Extract emissivity coefficients
        double j_s[4] = {};
        j_s[0] = j_i(m,n);
        j_s[1] = j_q(m,n);
        j_s[3] = j_v(m,n);

        // Extract absorptivity coefficients
        double alpha_s[4] = {};
        alpha_s[0] = alpha_i(m,n);
        alpha_s[1] = alpha_q(m,n);
        alpha_s[3] = alpha_v(m,n);

        // Extract rotativity coefficients
        double rho_s[4] = {};
        rho_s[1] = rho_q(m,n);
        rho_s[3] = rho_v(m,n);

        // Calculate optical depth
        double delta_tau = alpha_s[0] * delta_lambda_cgs;
        bool optically_thin = delta_tau <= delta_tau_max;

        // Prepare to couple to matter
        double alpha_sq = alpha_s[1] * alpha_s[1] + alpha_s[3] * alpha_s[3];
        double alpha_p = std::sqrt(alpha_sq);
        double rho_sq = rho_s[1] * rho_s[1] + rho_s[3] * rho_s[3];
        double rho_p = std::sqrt(rho_sq);
        double ss_end[4] = {};

        // Couple with no absorptivity or rotativity
        if (alpha_s[0] == 0.0 and rho_p == 0.0)
          for (int a = 0; a < 4; a++)
            ss_end[a] = ss_start[a] + j_s[a] * delta_lambda_cgs;

        // Couple with no polarized absorptivity or rotativity but with nonzero absorptivity
        else if (alpha_p == 0.0 and rho_p == 0.0)
        {
          // Optically thin case
          if (optically_thin)
          {
            double exp_neg = std::exp(-delta_tau);
            double expm1 = std::expm1(delta_tau);
            for (int a = 0; a < 4; a++)
              ss_end[a] = exp_neg * (ss_start[a] + j_s[a] / alpha_s[0] * expm1);
          }

          // Optically thick case
          else
            for (int a = 0; a < 4; a++)
              ss_end[a] = j_s[a] / alpha_s[0];
        }

        // Couple with no absorptivity but nonzero rotativity (I A2-A5)
        else if (alpha_s[0] == 0.0)
        {
          double cos_rho = std::cos(rho_p * delta_lambda_cgs);
          double sin_rho = std::sin(rho_p * delta_lambda_cgs);
          double sin_sq_rho = std::sin(rho_p * delta_lambda_cgs / 2.0);
          sin_sq_rho = sin_sq_rho * sin_sq_rho;
          double rho_ss = rho_s[1] * ss_start[1] + rho_s[3] * ss_start[3];
          ss_end[0] = ss_start[0];
          ss_end[1] = ss_start[1] * cos_rho + 2.0 * rho_s[1] * rho_ss / rho_sq * sin_sq_rho
              - rho_s[3] * ss_start[2] * sin_rho;
          ss_end[2] =
              ss_start[1] * cos_rho + (rho_s[3] * ss_start[1] - rho_s[1] * ss_start[3]) * sin_rho;
          ss_end[3] = ss_start[3] * cos_rho + 2.0 * rho_s[3] * rho_ss / rho_sq * sin_sq_rho
              + rho_s[1] * ss_start[2] * sin_rho;
          for (int a = 0; a < 4; a++)
            ss_end[a] += j_s[a] * delta_lambda_cgs;
        }

        // Couple with no rotativity but nonzero polarized absorptivity
        else if (rho_p == 0.0)
        {
          // Optically thin case (I A14-A17)
          if (optically_thin)
          {
            double exp_neg_i = std::exp(-delta_tau);
            double exp_neg_p = std::exp(-alpha_p * delta_lambda_cgs);
            double sinh_p = std::sinh(alpha_p * delta_lambda_cgs);
            double cosh_p = std::cosh(alpha_p * delta_lambda_cgs);
            double coshm1_p = 0.5 * (std::expm1(alpha_p * delta_lambda_cgs) + exp_neg_p - 1.0);
            double alpha_ss = alpha_s[1] * ss_start[1] + alpha_s[3] * ss_start[3];
            double alpha_j = alpha_s[1] * j_s[1] + alpha_s[3] * j_s[3];
            double alpha_i_p_factor = 1.0 / (alpha_s[0] * alpha_s[0] - alpha_sq);
            ss_end[0] = (ss_start[0] * cosh_p - alpha_ss / alpha_p * sinh_p) * exp_neg_i + alpha_j
                * alpha_i_p_factor * (-1.0 + (alpha_s[0] * sinh_p + alpha_p * cosh_p) / alpha_p
                * exp_neg_p) + alpha_s[0] * j_s[0] * alpha_i_p_factor * (1.0 - (alpha_s[0] * cosh_p
                + alpha_p * sinh_p) / alpha_s[0] * exp_neg_p);
            for (int a = 1; a < 4; a++)
            {
              double term_1 = (ss_start[a] + alpha_s[a] * alpha_ss / alpha_sq * coshm1_p
                  - ss_start[0] * alpha_s[a] / alpha_p * sinh_p) * exp_neg_i;
              double term_2 = j_s[a] * (1.0 - exp_neg_i) / alpha_s[0];
              double term_3 = alpha_j * alpha_s[a] / alpha_s[0] * alpha_i_p_factor * (1.0 - (1.0
                  - alpha_s[0] * alpha_s[0] / alpha_sq - alpha_s[0] / alpha_sq * (alpha_s[0]
                  * cosh_p + alpha_p * sinh_p)) * exp_neg_i);
              double term_4 = j_s[0] * alpha_s[a] / alpha_p * alpha_i_p_factor * (-alpha_p +
                  (alpha_p * cosh_p + alpha_s[0] * sinh_p) * exp_neg_i);
              ss_end[a] = term_1 + term_2 + term_3 + term_4;
            }
          }

          // Optically thick case
          else
          {
            double alpha_j = alpha_s[1] * j_s[1] + alpha_s[3] * j_s[3];
            ss_end[0] = (alpha_s[0] * j_s[0] - alpha_j) / (alpha_s[0] * alpha_s[0] - alpha_sq);
            for (int a = 1; a < 4; a++)
              ss_end[a] = (j_s[a] - alpha_s[a] * ss_end[0]) / alpha_s[0];
          }
        }

        // Couple with nonzero absorptivity and rotativity
        else
        {
          // Calculate coefficients needed for coupling matrices
          double alpha_rho = alpha_s[1] * rho_s[1] + alpha_s[3] * rho_s[3];
          double alpha_sq_rho_sq = alpha_sq - rho_sq;
          double lambda_a =
              std::sqrt(alpha_sq_rho_sq * alpha_sq_rho_sq / 4.0 + alpha_rho * alpha_rho);
          double lambda_b = alpha_sq_rho_sq / 2.0;
          double lambda_1 = std::sqrt(lambda_a + lambda_b);
          double lambda_2 = std::sqrt(lambda_a - lambda_b);
          double coefficient_theta = lambda_1 * lambda_1 + lambda_2 * lambda_2;
          double s = alpha_rho >= 0.0 ? 1.0 : -1.0;

          // Calculate coupling matrix 1
          double mm_1[4][4] = {};
          for (int a = 0; a < 4; a++)
            mm_1[a][a] = 1.0;

          // Calculate coupling matrix 2
          double mm_2[4][4] = {};
          mm_2[0][1] = lambda_2 * alpha_s[1] - s * lambda_1 * rho_s[1];
          mm_2[0][3] = lambda_2 * alpha_s[3] - s * lambda_1 * rho_s[3];
          mm_2[1][2] = s * lambda_1 * alpha_s[3] + lambda_2 * rho_s[3];
          mm_2[1][2] = s * lambda_1 * alpha_s[1] + lambda_2 * rho_s[1];
          mm_2[1][0] = mm_2[0][1];
          mm_2[2][0] = mm_2[0][2];
          mm_2[3][0] = mm_2[0][3];
          mm_2[2][1] = -mm_2[1][2];
          mm_2[3][1] = -mm_2[1][3];
          mm_2[3][2] = -mm_2[2][3];
          for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
              mm_2[a][b] *= 1.0 / coefficient_theta;

          // Calculate coupling matrix 3
          double mm_3[4][4] = {};
          mm_3[0][1] = lambda_1 * alpha_s[1] + s * lambda_2 * rho_s[1];
          mm_3[0][3] = lambda_1 * alpha_s[3] + s * lambda_2 * rho_s[3];
          mm_3[1][2] = -(s * lambda_2 * alpha_s[3] - lambda_1 * rho_s[3]);
          mm_3[1][2] = -(s * lambda_2 * alpha_s[1] - lambda_1 * rho_s[1]);
          mm_3[1][0] = mm_3[0][1];
          mm_3[2][0] = mm_3[0][2];
          mm_3[3][0] = mm_3[0][3];
          mm_3[2][1] = -mm_3[1][2];
          mm_3[3][1] = -mm_3[1][3];
          mm_3[3][2] = -mm_3[2][3];
          for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
              mm_3[a][b] *= 1.0 / coefficient_theta;

          // Calculate coupling matrix 4
          double mm_4[4][4] = {};
          mm_4[0][0] = (alpha_sq + rho_sq) / 2.0;
          mm_4[1][1] = alpha_s[1] * alpha_s[1] + rho_s[1] * rho_s[1] - (alpha_sq + rho_sq) / 2.0;
          mm_4[2][2] = -(alpha_sq + rho_sq) / 2.0;
          mm_4[3][3] = alpha_s[3] * alpha_s[3] + rho_s[3] * rho_s[3] - (alpha_sq + rho_sq) / 2.0;
          mm_4[0][2] = alpha_s[1] * rho_s[3] - alpha_s[3] * rho_s[1];
          mm_4[1][3] = alpha_s[3] * alpha_s[1] + rho_s[3] * rho_s[1];
          mm_4[1][0] = -mm_4[0][1];
          mm_4[2][0] = -mm_4[0][2];
          mm_4[3][0] = -mm_4[0][3];
          mm_4[2][1] = mm_4[1][2];
          mm_4[3][1] = mm_4[1][3];
          mm_4[3][2] = mm_4[2][3];
          for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
              mm_4[a][b] *= 2.0 / coefficient_theta;

          // Calculate coupling polynomial O (L 10)
          double exp, sin, cos, sinh, cosh;
          double oo[4][4] = {};
          if (optically_thin)
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
          double pp[4][4] = {};
          double f_1 = 1.0 / (alpha_s[0] * alpha_s[0] - lambda_1 * lambda_1);
          double f_2 = 1.0 / (alpha_s[0] * alpha_s[0] + lambda_2 * lambda_2);
          for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
            {
              double cosh_term =
                  -lambda_1 * f_1 * mm_3[a][b] + 0.5 * alpha_s[0] * f_1 * (mm_1[a][b] + mm_4[a][b]);
              double cos_term =
                  -lambda_2 * f_2 * mm_2[a][b] + 0.5 * alpha_s[0] * f_2 * (mm_1[a][b] - mm_4[a][b]);
              pp[a][b] = cosh_term + cos_term;
              if (optically_thin)
              {
                double sin_term = -alpha_s[0] * f_2 * mm_2[a][b]
                    - 0.5 * lambda_2 * f_2 * (mm_1[a][b] - mm_4[a][b]);
                double sinh_term = -alpha_s[0] * f_1 * mm_3[a][b]
                    + 0.5 * lambda_1 * f_1 * (mm_1[a][b] + mm_4[a][b]);
                pp[a][b] -=
                    exp * (cosh_term * cosh + cos_term * cos + sin_term * sin + sinh_term * sinh);
              }
            }


          // Apply coupling polynomials
          if (optically_thin)
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                ss_end[a] += pp[a][b] * j_s[b] + oo[a][b] * ss_start[b];
          else
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                ss_end[a] += pp[a][b] * j_s[b];
        }

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
            nn_con(m,mu,nu) = 0.0;
            for (int a = 0; a < 4; a++)
              for (int b = 0; b < 4; b++)
                nn_con(m,mu,nu) += tetrad(a,mu) * tetrad(b,nu) * nn_tet_con(a,b);
          }

        // Parallel-transport N by second half step
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
            nn_con_temp(mu,nu) = nn_con(m,mu,nu);
        for (int mu = 0; mu < 4; mu++)
          for (int nu = 0; nu < 4; nu++)
          {
            std::complex<double> dnn_dlambda = 0.0;
            for (int alpha = 0; alpha < 4; alpha++)
              for (int beta = 0; beta < 4; beta++)
                dnn_dlambda -= kcon[alpha] * (connection(mu,alpha,beta) * nn_con_temp(beta,nu)
                    + connection(nu,alpha,beta) * nn_con_temp(mu,beta));
            nn_con(m,mu,nu) += dnn_dlambda * (delta_lambda + delta_lambda_new) / 4.0;
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
    for (int m = 0; m < camera_num_pix; m++)
    {
      // Extract geodesic position and covariant momentum
      double x = camera_pos(m,1);
      double y = camera_pos(m,2);
      double z = camera_pos(m,3);
      double kcov[4];
      kcov[0] = camera_dir(m,0);
      kcov[1] = camera_dir(m,1);
      kcov[2] = camera_dir(m,2);
      kcov[3] = camera_dir(m,3);

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
                      * nn_con(m,alpha,beta);

      // Calculate orthonormal-frame Stokes quantities at camera location (I 14)
      image(0,m) = 0.5 * (nn_tet_cov(1,1) + nn_tet_cov(2,2)).real();
      image(1,m) = 0.5 * (nn_tet_cov(1,1) - nn_tet_cov(2,2)).real();
      image(2,m) = 0.5 * (nn_tet_cov(1,2) + nn_tet_cov(2,1)).real();
      image(3,m) = 0.5 * (nn_tet_cov(2,1) - nn_tet_cov(1,2)).imag();
    }

    // Transform invariant Stokes quantities (e.g. I_nu/nu^3) to standard ones (e.g. I_nu)
    double nu_cu = image_frequency * image_frequency * image_frequency;
    #pragma omp for schedule(static)
    for (int a = 0; a < 4; a++)
      for (int m = 0; m < camera_num_pix; m++)
        image(a,m) *= nu_cu;
  }
  return;
}
