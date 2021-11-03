// Blacklight geodesic integrator - camera definition

// C++ headers
#include <cmath>  // cos, exp, hypot, log, sin, sqrt

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "geodesic_integrator.hpp"
#include "../blacklight.hpp"        // enums
#include "../utils/array.hpp"       // Array

//--------------------------------------------------------------------------------------------------

// Function for setting up camera pixels and initial ray directions
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Allocates and initializes camera_pos[0], camera_dir[0], image_frequencies, and
//       momentum_factors[0].
//   Neglects spacetime curvature at camera location.
//   Symbols:
//     n: unit outward normal
//     u: unit right vector
//     v: unit up vector
void GeodesicIntegrator::InitializeCamera()
{
  // Calculate ideal image frequencies
  image_frequencies.Allocate(image_num_frequencies);
  if (image_num_frequencies == 1)
    image_frequencies(0) = image_frequency;
  else
  {
    image_frequencies(0) = image_frequency_start;
    image_frequencies(image_num_frequencies-1) = image_frequency_end;
    for (int l = 1; l < image_num_frequencies - 1; l++)
    {
      double frac = static_cast<double>(l) / static_cast<double>(image_num_frequencies - 1);
      if (image_frequency_spacing == FrequencySpacing::lin_freq)
        image_frequencies(l) =
            image_frequency_start + frac * (image_frequency_end - image_frequency_start);
      else if (image_frequency_spacing == FrequencySpacing::lin_wave)
        image_frequencies(l) = 1.0 / (1.0 / image_frequency_start
            + frac * (1.0 / image_frequency_end - 1.0 / image_frequency_start));
      else if (image_frequency_spacing == FrequencySpacing::log)
        image_frequencies(l) = std::exp(std::log(image_frequency_start)
            + frac * std::log(image_frequency_end / image_frequency_start));
    }
  }

  // Calculate trigonometric quantities
  double sth = std::sin(camera_th);
  double cth = std::cos(camera_th);
  double sph = std::sin(camera_ph);
  double cph = std::cos(camera_ph);
  double srot = std::sin(camera_rotation);
  double crot = std::cos(camera_rotation);

  // Calculate camera position
  cam_x[0] = 0.0;
  cam_x[1] = sth * (camera_r * cph - bh_a * sph);
  cam_x[2] = sth * (camera_r * sph + bh_a * cph);
  cam_x[3] = camera_r * cth;
  if (ray_flat)
  {
    cam_x[1] = camera_r * sth * cph;
    cam_x[2] = camera_r * sth * sph;
  }

  // Calculate metric in spherical coordinates
  double a2 = bh_a * bh_a;
  double r2 = camera_r * camera_r;
  double delta = r2 - 2.0 * bh_m * camera_r + a2;
  double sigma = r2 + a2 * cth * cth;
  double g_cov_r_r = 1.0 + 2.0 * bh_m * camera_r / sigma;
  double g_cov_r_th = 0.0;
  double g_cov_r_ph = -(1.0 + 2.0 * bh_m * camera_r / sigma) * bh_a * sth * sth;
  double g_cov_th_th = sigma;
  double g_cov_th_ph = 0.0;
  double g_cov_ph_ph = (r2 + a2 + 2.0 * bh_m * a2 * camera_r / sigma * sth * sth) * sth * sth;
  double g_con_t_t = -(1.0 + 2.0 * bh_m * camera_r / sigma);
  double g_con_t_r = 2.0 * bh_m * camera_r / sigma;
  double g_con_t_th = 0.0;
  double g_con_t_ph = 0.0;
  double g_con_r_r = delta / sigma;
  double g_con_r_th = 0.0;
  double g_con_r_ph = bh_a / sigma;
  double g_con_th_th = 1.0 / sigma;
  double g_con_th_ph = 0.0;
  double g_con_ph_ph = 1.0 / (sigma * sth * sth);
  if (ray_flat)
  {
    g_cov_r_r = 1.0;
    g_cov_r_ph = 0.0;
    g_cov_th_th = r2;
    g_cov_ph_ph = r2 * sth * sth;
    g_con_t_t = -1.0;
    g_con_t_r = 0.0;
    g_con_r_r = 1.0;
    g_con_r_ph = 0.0;
    g_con_th_th = 1.0 / r2;
    g_con_ph_ph = 1.0 / (r2 * sth * sth);
  }

  // Calculate camera velocity in spherical coordinates
  double alpha = 1.0 / std::sqrt(-g_con_t_t);
  double beta_con_r = -g_con_t_r / g_con_t_t;
  double beta_con_th = -g_con_t_th / g_con_t_t;
  double beta_con_ph = -g_con_t_ph / g_con_t_t;
  double utn = std::sqrt(1.0 + g_cov_r_r * camera_urn * camera_urn
      + 2.0 * g_cov_r_th * camera_urn * camera_uthn + 2.0 * g_cov_r_ph * camera_urn * camera_uphn
      + g_cov_th_th * camera_uthn * camera_uthn + 2.0 * g_cov_th_ph * camera_uthn * camera_uphn
      + g_cov_ph_ph * camera_uphn * camera_uphn);
  u_con[0] = utn / alpha;
  double ur = camera_urn - beta_con_r / alpha * utn;
  double uth = camera_uthn - beta_con_th / alpha * utn;
  double uph = camera_uphn - beta_con_ph / alpha * utn;

  // Calculate Jacobian of transformation
  double dx_dr = sth * cph;
  double dy_dr = sth * sph;
  double dz_dr = cth;
  double dx_dth = cth * (camera_r * cph - bh_a * sph);
  double dy_dth = cth * (camera_r * sph + bh_a * cph);
  double dz_dth = -camera_r * sth;
  double dx_dph = sth * (-camera_r * sph - bh_a * cph);
  double dy_dph = sth * (camera_r * cph - bh_a * sph);
  double dz_dph = 0.0;
  if (ray_flat)
  {
    dx_dth = camera_r * cth * cph;
    dy_dth = camera_r * cth * sph;
    dx_dph = -camera_r * sth * sph;
    dy_dph = camera_r * sth * cph;
  }

  // Calculate camera velocity
  u_con[1] = dx_dr * ur + dx_dth * uth + dx_dph * uph;
  u_con[2] = dy_dr * ur + dy_dth * uth + dy_dph * uph;
  u_con[3] = dz_dr * ur + dz_dth * uth + dz_dph * uph;
  double g_cov[4][4];
  CovariantGeodesicMetric(cam_x[1], cam_x[2], cam_x[3], g_cov);
  for (int mu = 0; mu < 4; mu++)
  {
    u_cov[mu] = 0.0;
    for (int nu = 0; nu < 4; nu++)
      u_cov[mu] += g_cov[mu][nu] * u_con[nu];
  }

  // Calculate photon momentum in spherical coordinates
  double g_con_rn_rn = (g_con_t_t * g_con_r_r - g_con_t_r * g_con_t_r) / g_con_t_t;
  double g_con_rn_thn = (g_con_t_t * g_con_r_th - g_con_t_r * g_con_t_th) / g_con_t_t;
  double g_con_rn_phn = (g_con_t_t * g_con_r_ph - g_con_t_r * g_con_t_ph) / g_con_t_t;
  double g_con_thn_thn = (g_con_t_t * g_con_th_th - g_con_t_th * g_con_t_th) / g_con_t_t;
  double g_con_thn_phn = (g_con_t_t * g_con_th_ph - g_con_t_th * g_con_t_ph) / g_con_t_t;
  double g_con_phn_phn = (g_con_t_t * g_con_ph_ph - g_con_t_ph * g_con_t_ph) / g_con_t_t;
  double k_rn = camera_k_r;
  double k_thn = camera_k_th;
  double k_phn = camera_k_ph;
  double k_tn = -std::sqrt(g_con_rn_rn * k_rn * k_rn + 2.0 * g_con_rn_thn * k_rn * k_thn
      + 2.0 * g_con_rn_phn * k_rn * k_phn + g_con_thn_thn * k_thn * k_thn
      + 2.0 * g_con_thn_phn * k_thn * k_phn + g_con_phn_phn * k_phn * k_phn);
  double k_t = alpha * k_tn + (beta_con_r * k_rn + beta_con_th * k_thn + beta_con_ph * k_phn);

  // Calculate Jacobian of transformation
  double rr2 = cam_x[1] * cam_x[1] + cam_x[2] * cam_x[2] + cam_x[3] * cam_x[3];
  double dr_dx = camera_r * cam_x[1] / (2.0 * r2 - rr2 + a2);
  double dr_dy = camera_r * cam_x[2] / (2.0 * r2 - rr2 + a2);
  double dr_dz = (camera_r * cam_x[3] + a2 * cam_x[3] / camera_r) / (2.0 * r2 - rr2 + a2);
  double dth_dx = cam_x[3] * dr_dx / (r2 * sth);
  double dth_dy = cam_x[3] * dr_dy / (r2 * sth);
  double dth_dz = (cam_x[3] * dr_dz - camera_r) / (r2 * sth);
  double dph_dx =
      -cam_x[2] / (cam_x[1] * cam_x[1] + cam_x[2] * cam_x[2]) + bh_a / (r2 + a2) * dr_dx;
  double dph_dy = cam_x[1] / (cam_x[1] * cam_x[1] + cam_x[2] * cam_x[2]) + bh_a / (r2 + a2) * dr_dy;
  double dph_dz = bh_a / (r2 + a2) * dr_dz;
  if (ray_flat)
  {
    dr_dx = cam_x[1] / camera_r;
    dr_dy = cam_x[2] / camera_r;
    dr_dz = cam_x[3] / camera_r;
    dth_dx = cth * cph / camera_r;
    dth_dy = cth * sph / camera_r;
    dth_dz = -sth / camera_r;
    dph_dx = -sph / (camera_r * sth);
    dph_dy = cph / (camera_r * sth);
    dph_dz = 0.0;
  }

  // Calculate photon momentum
  double k_x = dr_dx * camera_k_r + dth_dx * camera_k_th + dph_dx * camera_k_ph;
  double k_y = dr_dy * camera_k_r + dth_dy * camera_k_th + dph_dy * camera_k_ph;
  double k_z = dr_dz * camera_k_r + dth_dz * camera_k_th + dph_dz * camera_k_ph;
  double k_tc = u_con[0] * k_t + u_con[1] * k_x + u_con[2] * k_y + u_con[3] * k_z;

  // Calculate contravariant metric in camera frame
  double g_con[4][4];
  ContravariantGeodesicMetric(cam_x[1], cam_x[2], cam_x[3], g_con);
  double g_con_xc_xc = g_con[1][1] + u_con[1] * u_con[1];
  double g_con_xc_yc = g_con[1][2] + u_con[1] * u_con[2];
  double g_con_xc_zc = g_con[1][3] + u_con[1] * u_con[3];
  double g_con_yc_yc = g_con[2][2] + u_con[2] * u_con[2];
  double g_con_yc_zc = g_con[2][3] + u_con[2] * u_con[3];
  double g_con_zc_zc = g_con[3][3] + u_con[3] * u_con[3];

  // Calculate camera normal direction in camera frame
  double norm_cov_xc = k_x - u_cov[1] / u_cov[0] * k_t;
  double norm_cov_yc = k_y - u_cov[2] / u_cov[0] * k_t;
  double norm_cov_zc = k_z - u_cov[3] / u_cov[0] * k_t;
  norm_con_c[0] = -k_tc;
  norm_con_c[1] = g_con_xc_xc * norm_cov_xc + g_con_xc_yc * norm_cov_yc + g_con_xc_zc * norm_cov_zc;
  norm_con_c[2] = g_con_xc_yc * norm_cov_xc + g_con_yc_yc * norm_cov_yc + g_con_yc_zc * norm_cov_zc;
  norm_con_c[3] = g_con_xc_zc * norm_cov_xc + g_con_yc_zc * norm_cov_yc + g_con_zc_zc * norm_cov_zc;
  double norm_norm = std::sqrt(norm_cov_xc * norm_con_c[1] + norm_cov_yc * norm_con_c[2]
      + norm_cov_zc * norm_con_c[3]);
  norm_cov_xc /= norm_norm;
  norm_cov_yc /= norm_norm;
  norm_cov_zc /= norm_norm;
  norm_con_c[0] /= norm_norm;
  norm_con_c[1] /= norm_norm;
  norm_con_c[2] /= norm_norm;
  norm_con_c[3] /= norm_norm;
  norm_con[0] = u_con[0] * norm_con_c[0]
      - (u_cov[1] * norm_con_c[1] + u_cov[2] * norm_con_c[2] + u_cov[3] * norm_con_c[3]) / u_cov[0];
  norm_con[1] = norm_con_c[1] + u_con[1] * norm_con_c[0];
  norm_con[2] = norm_con_c[2] + u_con[2] * norm_con_c[0];
  norm_con[3] = norm_con_c[3] + u_con[3] * norm_con_c[0];

  // Define unprojected vertical direction in camera frame
  double up_con_xc = 0.0;
  double up_con_yc = 0.0;
  double up_con_zc = 1.0;
  if (camera_pole)
  {
    up_con_yc = 1.0;
    up_con_zc = 0.0;
  }

  // Calculate covariant metric in camera frame
  double g_cov_xc_xc = g_cov[1][1] - u_cov[1] / u_cov[0] * g_cov[1][0]
      - u_cov[1] / u_cov[0] * g_cov[1][0]
      + u_cov[1] * u_cov[1] / (u_cov[0] * u_cov[0]) * g_cov[0][0];
  double g_cov_xc_yc = g_cov[1][2] - u_cov[1] / u_cov[0] * g_cov[2][0]
      - u_cov[2] / u_cov[0] * g_cov[1][0]
      + u_cov[1] * u_cov[2] / (u_cov[0] * u_cov[0]) * g_cov[0][0];
  double g_cov_xc_zc = g_cov[1][3] - u_cov[1] / u_cov[0] * g_cov[3][0]
      - u_cov[3] / u_cov[0] * g_cov[1][0]
      + u_cov[1] * u_cov[3] / (u_cov[0] * u_cov[0]) * g_cov[0][0];
  double g_cov_yc_yc = g_cov[2][2] - u_cov[2] / u_cov[0] * g_cov[2][0]
      - u_cov[2] / u_cov[0] * g_cov[2][0]
      + u_cov[2] * u_cov[2] / (u_cov[0] * u_cov[0]) * g_cov[0][0];
  double g_cov_yc_zc = g_cov[2][3] - u_cov[2] / u_cov[0] * g_cov[3][0]
      - u_cov[3] / u_cov[0] * g_cov[2][0]
      + u_cov[2] * u_cov[3] / (u_cov[0] * u_cov[0]) * g_cov[0][0];
  double g_cov_zc_zc = g_cov[3][3] - u_cov[3] / u_cov[0] * g_cov[3][0]
      - u_cov[3] / u_cov[0] * g_cov[3][0]
      + u_cov[3] * u_cov[3] / (u_cov[0] * u_cov[0]) * g_cov[0][0];

  // Calculate camera vertical direction without rotation in camera frame
  double up_norm = up_con_xc * norm_cov_xc + up_con_yc * norm_cov_yc + up_con_zc * norm_cov_zc;
  vert_con_c[0] = 0.0;
  vert_con_c[1] = up_con_xc - up_norm * norm_con_c[1];
  vert_con_c[2] = up_con_yc - up_norm * norm_con_c[2];
  vert_con_c[3] = up_con_zc - up_norm * norm_con_c[3];
  double vert_cov_xc =
      g_cov_xc_xc * vert_con_c[1] + g_cov_xc_yc * vert_con_c[2] + g_cov_xc_zc * vert_con_c[3];
  double vert_cov_yc =
      g_cov_xc_yc * vert_con_c[1] + g_cov_yc_yc * vert_con_c[2] + g_cov_yc_zc * vert_con_c[3];
  double vert_cov_zc =
      g_cov_xc_zc * vert_con_c[1] + g_cov_yc_zc * vert_con_c[2] + g_cov_zc_zc * vert_con_c[3];
  double vert_norm = std::sqrt(vert_cov_xc * vert_con_c[1] + vert_cov_yc * vert_con_c[2]
      + vert_cov_zc * vert_con_c[3]);
  vert_cov_xc /= vert_norm;
  vert_cov_yc /= vert_norm;
  vert_cov_zc /= vert_norm;
  vert_con_c[1] /= vert_norm;
  vert_con_c[2] /= vert_norm;
  vert_con_c[3] /= vert_norm;

  // Calculate determinant of metric in camera frame
  double det = g_cov_xc_xc * (g_cov_yc_yc * g_cov_zc_zc - g_cov_yc_zc * g_cov_yc_zc)
      + g_cov_xc_yc * (g_cov_yc_zc * g_cov_xc_zc - g_cov_xc_yc * g_cov_zc_zc)
      + g_cov_xc_zc * (g_cov_xc_yc * g_cov_yc_zc - g_cov_yc_yc * g_cov_xc_zc);
  double det_sqrt = std::sqrt(det);

  // Calculate camera horizontal direction without rotation in camera frame
  hor_con_c[0] = 0.0;
  hor_con_c[1] = (vert_cov_yc * norm_cov_zc - vert_cov_zc * norm_cov_yc) / det_sqrt;
  hor_con_c[2] = (vert_cov_zc * norm_cov_xc - vert_cov_xc * norm_cov_zc) / det_sqrt;
  hor_con_c[3] = (vert_cov_xc * norm_cov_yc - vert_cov_yc * norm_cov_xc) / det_sqrt;

  // Calculate camera direction with rotation in camera frame
  double temp_hor_con_xc = hor_con_c[1];
  double temp_hor_con_yc = hor_con_c[2];
  double temp_hor_con_zc = hor_con_c[3];
  double temp_vert_con_xc = vert_con_c[1];
  double temp_vert_con_yc = vert_con_c[2];
  double temp_vert_con_zc = vert_con_c[3];
  hor_con_c[1] = temp_hor_con_xc * crot - temp_vert_con_xc * srot;
  hor_con_c[2] = temp_hor_con_yc * crot - temp_vert_con_yc * srot;
  hor_con_c[3] = temp_hor_con_zc * crot - temp_vert_con_zc * srot;
  vert_con_c[1] = temp_vert_con_xc * crot + temp_hor_con_xc * srot;
  vert_con_c[2] = temp_vert_con_yc * crot + temp_hor_con_yc * srot;
  vert_con_c[3] = temp_vert_con_zc * crot + temp_hor_con_zc * srot;

  // Allocate arrays
  camera_pos[0].Allocate(camera_num_pix, 4);
  camera_dir[0].Allocate(camera_num_pix, 4);
  momentum_factors[0].Allocate(camera_num_pix);

  // Initialize plane-parallel camera
  if (camera_type == Camera::plane)
  {
    #pragma omp parallel for schedule(static)
    for (int m = 0; m < camera_num_pix; m++)
    {
      int m2 = m / camera_resolution;
      int m1 = m % camera_resolution;
      double u_ind = (m1 - camera_resolution / 2.0 + 0.5) / camera_resolution;
      double v_ind = (m2 - camera_resolution / 2.0 + 0.5) / camera_resolution;
      SetPixelPlane(u_ind, v_ind, m, camera_pos[0], camera_dir[0], momentum_factors[0]);
    }
  }

  // Initialize pinhole camera
  if (camera_type == Camera::pinhole)
  {
    #pragma omp parallel for schedule(static)
    for (int m = 0; m < camera_num_pix; m++)
    {
      int m2 = m / camera_resolution;
      int m1 = m % camera_resolution;
      double u_ind = (m1 - camera_resolution / 2.0 + 0.5) / camera_resolution;
      double v_ind = (m2 - camera_resolution / 2.0 + 0.5) / camera_resolution;
      SetPixelPinhole(u_ind, v_ind, m, camera_pos[0], camera_dir[0], momentum_factors[0]);
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for adaptively adding more pixels to camera
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes block_count[adaptive_level] and block_count[adaptive_level-1] have been set.
//   Allocates and initializes camera_loc[adaptive_level], camera_pos[adaptive_level],
//       camera_dir[adaptive_level], and momentum_factors[adaptive_level].
void GeodesicIntegrator::AugmentCamera()
{
  // Allocate storage for new blocks
  int block_count = block_counts[adaptive_level];
  camera_loc[adaptive_level].Allocate(block_count, 2);
  camera_pos[adaptive_level].Allocate(block_count * block_num_pix, 4);
  camera_dir[adaptive_level].Allocate(block_count * block_num_pix, 4);
  momentum_factors[adaptive_level].Allocate(block_count * block_num_pix);

  // Prepare to go through blocks
  int block_count_old = block_counts[adaptive_level-1];
  int effective_resolution = camera_resolution;
  for (int n = 1; n <= adaptive_level; n++)
    effective_resolution *= 2;
  Array<double> camera_pos_block;
  Array<double> camera_dir_block;
  Array<double> momentum_factors_block;

  // Go through blocks at previous level
  for (int block_old = 0, block = 0; block_old < block_count_old; block_old++)
    if (refinement_flags[adaptive_level-1](block_old))
    {
      // Locate block
      int block_v_old = camera_loc[adaptive_level-1](block_old,0);
      int block_u_old = camera_loc[adaptive_level-1](block_old,1);

      // Go through new blocks at current level
      for (int block_v = 2 * block_v_old; block_v <= 2 * block_v_old + 1; block_v++)
        for (int block_u = 2 * block_u_old; block_u <= 2 * block_u_old + 1; block_u++, block++)
        {
          // Calculate location in image plane
          camera_loc[adaptive_level](block,0) = block_v;
          camera_loc[adaptive_level](block,1) = block_u;
          camera_pos_block = camera_pos[adaptive_level];
          camera_dir_block = camera_dir[adaptive_level];
          momentum_factors_block = momentum_factors[adaptive_level];
          camera_pos_block.Slice(2, block * block_num_pix, block * block_num_pix);
          camera_dir_block.Slice(2, block * block_num_pix, block * block_num_pix);
          momentum_factors_block.Slice(1, block * block_num_pix, block * block_num_pix);
          int m_offset = block_v * adaptive_block_size;
          int l_offset = block_u * adaptive_block_size;

          // Initialize position and direction for plane-parallel camera
          if (camera_type == Camera::plane)
          {
            #pragma omp parallel for schedule(static)
            for (int m = 0; m < block_num_pix; m++)
            {
              int m2 = m / adaptive_block_size;
              int m1 = m % adaptive_block_size;
              double u_ind =
                  (m1 + l_offset - effective_resolution / 2.0 + 0.5) / effective_resolution;
              double v_ind =
                  (m2 + m_offset - effective_resolution / 2.0 + 0.5) / effective_resolution;
              SetPixelPlane(u_ind, v_ind, m, camera_pos_block, camera_dir_block,
                  momentum_factors_block);
            }
          }

          // Initialize position and direction for pinhole camera
          if (camera_type == Camera::pinhole)
          {
            #pragma omp parallel for schedule(static)
            for (int m = 0; m < block_num_pix; m++)
            {
              int m2 = m / adaptive_block_size;
              int m1 = m % adaptive_block_size;
              double u_ind =
                  (m1 + l_offset - effective_resolution / 2.0 + 0.5) / effective_resolution;
              double v_ind =
                  (m2 + m_offset - effective_resolution / 2.0 + 0.5) / effective_resolution;
              SetPixelPinhole(u_ind, v_ind, m, camera_pos_block, camera_dir_block,
                  momentum_factors_block);
            }
          }
        }
    }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating individual pixel position and direction for plane-parallel camera
// Inputs:
//   u_ind: fractional horizontal coordinate, between -0.5 (left edge) and +0.5 (right edge)
//   v_ind: fractional vertical coordinate, between -0.5 (bottom edge) and +0.5 (top edge)
//   m: index of pixel in arrays, corresponding to the second-to-last dimension
// Outputs:
//   position: appropriate values of array updated with spacetime location of pixel
//   direction: appropriate values of array updated with contravariant spatial momentum of light
//       seen by pixel
//   factor: appropriate values of array updated, as function of pixel and ideal frequency (at
//       preferred normalization location), with factor that, when multiplied by momentum,
//       appropriately normalizes the ray
// Notes:
//   Assumes cam_x, u_con, u_cov, norm_con, hor_con_c, vert_con_c, and image_frequencies have been
//       set.
//   Quadratic solved as follows:
//     Outside ergosphere: unique positive root.
//     On ergosphere, assuming g_{0i} p^i < 0: unique root, which will be positive
//     Inside ergosphere, assuming g_{0i} p^i < 0: lesser positive root, which remains finite as
//         ergosphere is approached.
void GeodesicIntegrator::SetPixelPlane(double u_ind, double v_ind, int m, Array<double> &position,
    Array<double> &direction, Array<double> &factor)
{
  // Set pixel position
  double u = u_ind * bh_m * camera_width;
  double v = v_ind * bh_m * camera_width;
  double dtc = u * hor_con_c[0] + v * vert_con_c[0];
  double dxc = u * hor_con_c[1] + v * vert_con_c[1];
  double dyc = u * hor_con_c[2] + v * vert_con_c[2];
  double dzc = u * hor_con_c[3] + v * vert_con_c[3];
  double dt = u_con[0] * dtc - (u_cov[1] * dxc + u_cov[2] * dyc + u_cov[3] * dzc) / u_cov[0];
  double dx = dxc + u_con[1] * dtc;
  double dy = dyc + u_con[2] * dtc;
  double dz = dzc + u_con[3] * dtc;
  position(m,0) = cam_x[0] + dt;
  position(m,1) = cam_x[1] + dx;
  position(m,2) = cam_x[2] + dy;
  position(m,3) = cam_x[3] + dz;

  // Calculate pixel direction
  double p[4];
  p[1] = norm_con[1];
  p[2] = norm_con[2];
  p[3] = norm_con[3];

  // Calculate time component of momentum
  double gcov[4][4];
  CovariantGeodesicMetric(position(m,1), position(m,2), position(m,3), gcov);
  double temp_a = gcov[0][0];
  double temp_b = 0.0;
  for (int a = 1; a < 4; a++)
    temp_b += 2.0 * gcov[0][a] * p[a];
  double temp_c = 0.0;
  for (int a = 1; a < 4; a++)
    for (int b = 1; b < 4; b++)
      temp_c += gcov[a][b] * p[a] * p[b];
  double temp_d = std::sqrt(std::max(temp_b * temp_b - 4.0 * temp_a * temp_c, 0.0));
  p[0] = temp_a == 0.0 ? -temp_c / (2.0 * temp_b)
      : (temp_b < 0.0 ? 2.0 * temp_c / (temp_d - temp_b) : -(temp_b + temp_d) / (2.0 * temp_a));

  // Set pixel momentum
  for (int mu = 0; mu < 4; mu++)
  {
    direction(m,mu) = 0.0;
    for (int nu = 0; nu < 4; nu++)
      direction(m,mu) += gcov[mu][nu] * p[nu];
  }

  // Set momentum factor
  double nu_local = 0.0;
  if (image_normalization == FrequencyNormalization::camera)
    for (int mu = 0; mu < 4; mu++)
      nu_local -= direction(m,mu) * u_con[mu];
  else if (image_normalization == FrequencyNormalization::infinity)
    nu_local = -direction(m,0);
  factor(m) = 1.0 / nu_local;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating individual pixel position and direction for pinhole camera
// Inputs:
//   u_ind: fractional horizontal coordinate, between -0.5 (left edge) and +0.5 (right edge)
//   v_ind: fractional vertical coordinate, between -0.5 (bottom edge) and +0.5 (top edge)
//   m: index of pixel in arrays, corresponding to the second-to-last dimension
// Outputs:
//   position: appropriate values of array updated with spacetime location of pixel
//   direction: appropriate values of array updated with contravariant spatial momentum of light
//       seen by pixel
//   factor: appropriate values of array updated, as function of pixel and ideal frequency (at
//       preferred normalization location), with factor that, when multiplied by momentum,
//       appropriately normalizes the ray
// Notes:
//   Assumes cam_x, u_con, norm_con_c, hor_con_c, vert_con_c, and image_frequencies have been set.
//   Quadratic solved as follows:
//     Outside ergosphere: unique positive root.
//     On ergosphere, assuming g_{0i} p^i < 0: unique root, which will be positive
//     Inside ergosphere, assuming g_{0i} p^i < 0: lesser positive root, which remains finite as
//         ergosphere is approached.
void GeodesicIntegrator::SetPixelPinhole(double u_ind, double v_ind, int m, Array<double> &position,
    Array<double> &direction, Array<double> &factor)
{
  // Set pixel position
  position(m,0) = cam_x[0];
  position(m,1) = cam_x[1];
  position(m,2) = cam_x[2];
  position(m,3) = cam_x[3];

  // Calculate pixel direction
  double u = u_ind * bh_m * camera_width;
  double v = v_ind * bh_m * camera_width;
  double normalization = std::hypot(u, v, camera_r);
  double frac_norm = camera_r / normalization;
  double frac_hor = -u / normalization;
  double frac_vert = -v / normalization;
  double dir_con_tc = norm_con_c[0];
  double dir_con_xc =
      frac_norm * norm_con_c[1] + frac_hor * hor_con_c[1] + frac_vert * vert_con_c[1];
  double dir_con_yc =
      frac_norm * norm_con_c[2] + frac_hor * hor_con_c[2] + frac_vert * vert_con_c[2];
  double dir_con_zc =
      frac_norm * norm_con_c[3] + frac_hor * hor_con_c[3] + frac_vert * vert_con_c[3];
  double dir_con_x = dir_con_xc + u_con[1] * dir_con_tc;
  double dir_con_y = dir_con_yc + u_con[2] * dir_con_tc;
  double dir_con_z = dir_con_zc + u_con[3] * dir_con_tc;
  double p[4];
  p[1] = dir_con_x;
  p[2] = dir_con_y;
  p[3] = dir_con_z;

  // Calculate time component of momentum
  double gcov[4][4];
  CovariantGeodesicMetric(position(m,1), position(m,2), position(m,3), gcov);
  double temp_a = gcov[0][0];
  double temp_b = 0.0;
  for (int a = 1; a < 4; a++)
    temp_b += 2.0 * gcov[0][a] * p[a];
  double temp_c = 0.0;
  for (int a = 1; a < 4; a++)
    for (int b = 1; b < 4; b++)
      temp_c += gcov[a][b] * p[a] * p[b];
  double temp_d = std::sqrt(std::max(temp_b * temp_b - 4.0 * temp_a * temp_c, 0.0));
  p[0] = temp_a == 0.0 ? -temp_c / (2.0 * temp_b)
      : (temp_b < 0.0 ? 2.0 * temp_c / (temp_d - temp_b) : -(temp_b + temp_d) / (2.0 * temp_a));

  // Set pixel momentum
  for (int mu = 0; mu < 4; mu++)
  {
    direction(m,mu) = 0.0;
    for (int nu = 0; nu < 4; nu++)
      direction(m,mu) += gcov[mu][nu] * p[nu];
  }

  // Set momentum factor
  double nu_local = 0.0;
  if (image_normalization == FrequencyNormalization::camera)
    for (int mu = 0; mu < 4; mu++)
      nu_local -= direction(m,mu) * u_con[mu];
  else if (image_normalization == FrequencyNormalization::infinity)
    nu_local = -direction(m,0);
  factor(m) = 1.0 / nu_local;
  return;
}
