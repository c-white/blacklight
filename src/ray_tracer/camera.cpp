// Blacklight ray tracer - camera definition

// C++ headers
#include <cmath>  // cos, hypot, sin, sqrt

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "ray_tracer.hpp"
#include "../blacklight.hpp"   // enums
#include "../utils/array.hpp"  // Array

//--------------------------------------------------------------------------------------------------

// Function for setting up camera pixels and initial ray directions
// Inputs: (none)
// Output: (none)
// Notes:
//   Allocates and initializes image_position and image_direction except for time components of
//       image_direction.
//   Neglects spacetime curvature at camera location.
//   Symbols:
//     n: unit outward normal
//     u: unit right vector
//     v: unit up vector
void RayTracer::InitializeCamera()
{
  // Calculate trigonometric quantities
  double sth = std::sin(image_th);
  double cth = std::cos(image_th);
  double sph = std::sin(image_ph);
  double cph = std::cos(image_ph);
  double srot = std::sin(image_rotation);
  double crot = std::cos(image_rotation);

  // Calculate camera position
  double t = 0.0;
  double x = sth * (image_r * cph - bh_a * sph);
  double y = sth * (image_r * sph + bh_a * cph);
  double z = image_r * cth;
  if (ray_flat)
  {
    x = image_r * sth * cph;
    y = image_r * sth * sph;
  }

  // Calculate metric in spherical coordinates
  double a2 = bh_a * bh_a;
  double r2 = image_r * image_r;
  double delta = r2 - 2.0 * bh_m * image_r + a2;
  double sigma = r2 + a2 * cth * cth;
  double gcov_r_r = 1.0 + 2.0 * bh_m * image_r / sigma;
  double gcov_r_th = 0.0;
  double gcov_r_ph = -(1.0 + 2.0 * bh_m * image_r / sigma) * bh_a * sth * sth;
  double gcov_th_th = sigma;
  double gcov_th_ph = 0.0;
  double gcov_ph_ph = (r2 + a2 + 2.0 * bh_m * a2 * image_r / sigma * sth * sth) * sth * sth;
  double gcon_t_t = -(1.0 + 2.0 * bh_m * image_r / sigma);
  double gcon_t_r = 2.0 * bh_m * image_r / sigma;
  double gcon_t_th = 0.0;
  double gcon_t_ph = 0.0;
  double gcon_r_r = delta / sigma;
  double gcon_r_th = 0.0;
  double gcon_r_ph = bh_a / sigma;
  double gcon_th_th = 1.0 / sigma;
  double gcon_th_ph = 0.0;
  double gcon_ph_ph = 1.0 / (sigma * sth * sth);
  if (ray_flat)
  {
    gcov_r_r = 1.0;
    gcov_r_ph = 0.0;
    gcov_th_th = r2;
    gcov_ph_ph = r2 * sth * sth;
    gcon_t_t = -1.0;
    gcon_t_r = 0.0;
    gcon_r_r = 1.0;
    gcon_r_ph = 0.0;
    gcon_th_th = 1.0 / r2;
    gcon_ph_ph = 1.0 / (r2 * sth * sth);
  }

  // Calculate camera velocity in spherical coordinates
  double alpha = 1.0 / std::sqrt(-gcon_t_t);
  double beta_con_r = -gcon_t_r / gcon_t_t;
  double beta_con_th = -gcon_t_th / gcon_t_t;
  double beta_con_ph = -gcon_t_ph / gcon_t_t;
  double utn = std::sqrt(1.0 + gcov_r_r * image_urn * image_urn
      + 2.0 * gcov_r_th * image_urn * image_uthn + 2.0 * gcov_r_ph * image_urn * image_uphn
      + gcov_th_th * image_uthn * image_uthn + 2.0 * gcov_th_ph * image_uthn * image_uphn
      + gcov_ph_ph * image_uphn * image_uphn);
  double ut = utn / alpha;
  double ur = image_urn - beta_con_r / alpha * utn;
  double uth = image_uthn - beta_con_th / alpha * utn;
  double uph = image_uphn - beta_con_ph / alpha * utn;

  // Calculate Jacobian of transformation
  double dx_dr = sth * cph;
  double dy_dr = sth * sph;
  double dz_dr = cth;
  double dx_dth = cth * (image_r * cph - bh_a * sph);
  double dy_dth = cth * (image_r * sph + bh_a * cph);
  double dz_dth = -image_r * sth;
  double dx_dph = sth * (-image_r * sph - bh_a * cph);
  double dy_dph = sth * (image_r * cph - bh_a * sph);
  double dz_dph = 0.0;
  if (ray_flat)
  {
    dx_dth = image_r * cth * cph;
    dy_dth = image_r * cth * sph;
    dx_dph = -image_r * sth * sph;
    dy_dph = image_r * sth * cph;
  }

  // Calculate camera velocity
  double ux = dx_dr * ur + dx_dth * uth + dx_dph * uph;
  double uy = dy_dr * ur + dy_dth * uth + dy_dph * uph;
  double uz = dz_dr * ur + dz_dth * uth + dz_dph * uph;
  Array<double> gcov(4, 4);
  CovariantGeodesicMetric(x, y, z, gcov);
  double u_t = gcov(0,0) * ut + gcov(0,1) * ux + gcov(0,2) * uy + gcov(0,3) * uz;
  double u_x = gcov(1,0) * ut + gcov(1,1) * ux + gcov(1,2) * uy + gcov(1,3) * uz;
  double u_y = gcov(2,0) * ut + gcov(2,1) * ux + gcov(2,2) * uy + gcov(2,3) * uz;
  double u_z = gcov(3,0) * ut + gcov(3,1) * ux + gcov(3,2) * uy + gcov(3,3) * uz;

  // Calculate photon momentum in spherical coordinates
  double gcon_rn_rn = (gcon_t_t * gcon_r_r - gcon_t_r * gcon_t_r) / gcon_t_t;
  double gcon_rn_thn = (gcon_t_t * gcon_r_th - gcon_t_r * gcon_t_th) / gcon_t_t;
  double gcon_rn_phn = (gcon_t_t * gcon_r_ph - gcon_t_r * gcon_t_ph) / gcon_t_t;
  double gcon_thn_thn = (gcon_t_t * gcon_th_th - gcon_t_th * gcon_t_th) / gcon_t_t;
  double gcon_thn_phn = (gcon_t_t * gcon_th_ph - gcon_t_th * gcon_t_ph) / gcon_t_t;
  double gcon_phn_phn = (gcon_t_t * gcon_ph_ph - gcon_t_ph * gcon_t_ph) / gcon_t_t;
  double k_rn = image_k_r;
  double k_thn = image_k_th;
  double k_phn = image_k_ph;
  double k_tn = -std::sqrt(gcon_rn_rn * k_rn * k_rn + 2.0 * gcon_rn_thn * k_rn * k_thn
      + 2.0 * gcon_rn_phn * k_rn * k_phn + gcon_thn_thn * k_thn * k_thn
      + 2.0 * gcon_thn_phn * k_thn * k_phn + gcon_phn_phn * k_phn * k_phn);
  double k_t = alpha * k_tn + (beta_con_r * k_rn + beta_con_th * k_thn + beta_con_ph * k_phn);

  // Calculate Jacobian of transformation
  double rr2 = x * x + y * y + z * z;
  double dr_dx = image_r * x / (2.0 * r2 - rr2 + a2);
  double dr_dy = image_r * y / (2.0 * r2 - rr2 + a2);
  double dr_dz = (image_r * z + a2 * z / image_r) / (2.0 * r2 - rr2 + a2);
  double dth_dx = z * dr_dx / (r2 * sth);
  double dth_dy = z * dr_dy / (r2 * sth);
  double dth_dz = (z * dr_dz - image_r) / (r2 * sth);
  double dph_dx = -y / (x * x + y * y) + bh_a / (r2 + a2) * dr_dx;
  double dph_dy = x / (x * x + y * y) + bh_a / (r2 + a2) * dr_dy;
  double dph_dz = bh_a / (r2 + a2) * dr_dz;
  if (ray_flat)
  {
    dr_dx = x / image_r;
    dr_dy = y / image_r;
    dr_dz = z / image_r;
    dth_dx = cth * cph / image_r;
    dth_dy = cth * sph / image_r;
    dth_dz = -sth / image_r;
    dph_dx = -sph / (image_r * sth);
    dph_dy = cph / (image_r * sth);
    dph_dz = 0.0;
  }

  // Calculate photon momentum
  double k_x = dr_dx * image_k_r + dth_dx * image_k_th + dph_dx * image_k_ph;
  double k_y = dr_dy * image_k_r + dth_dy * image_k_th + dph_dy * image_k_ph;
  double k_z = dr_dz * image_k_r + dth_dz * image_k_th + dph_dz * image_k_ph;
  double k_tc = ut * k_t + ux * k_x + uy * k_y + uz * k_z;

  // Calculate momentum normalization
  switch (image_normalization)
  {
    case FrequencyNormalization::camera:
    {
      momentum_factor = -image_frequency / k_tc;
      break;
    }
    case FrequencyNormalization::infinity:
    {
      momentum_factor = -image_frequency / k_t;
      break;
    }
  }

  // Calculate contravariant metric in camera frame
  Array<double> gcon(4, 4);
  ContravariantGeodesicMetric(x, y, z, gcon);
  double gcon_xc_xc = gcon(1,1) + ux * ux;
  double gcon_xc_yc = gcon(1,2) + ux * uy;
  double gcon_xc_zc = gcon(1,3) + ux * uz;
  double gcon_yc_yc = gcon(2,2) + uy * uy;
  double gcon_yc_zc = gcon(2,3) + uy * uz;
  double gcon_zc_zc = gcon(3,3) + uz * uz;

  // Calculate camera normal direction in camera frame
  double norm_cov_xc = k_x - u_x / u_t * k_t;
  double norm_cov_yc = k_y - u_y / u_t * k_t;
  double norm_cov_zc = k_z - u_z / u_t * k_t;
  double norm_con_tc = -k_tc;
  double norm_con_xc =
      gcon_xc_xc * norm_cov_xc + gcon_xc_yc * norm_cov_yc + gcon_xc_zc * norm_cov_zc;
  double norm_con_yc =
      gcon_xc_yc * norm_cov_xc + gcon_yc_yc * norm_cov_yc + gcon_yc_zc * norm_cov_zc;
  double norm_con_zc =
      gcon_xc_zc * norm_cov_xc + gcon_yc_zc * norm_cov_yc + gcon_zc_zc * norm_cov_zc;
  double norm_norm =
      std::sqrt(norm_cov_xc * norm_con_xc + norm_cov_yc * norm_con_yc + norm_cov_zc * norm_con_zc);
  norm_cov_xc /= norm_norm;
  norm_cov_yc /= norm_norm;
  norm_cov_zc /= norm_norm;
  norm_con_tc /= norm_norm;
  norm_con_xc /= norm_norm;
  norm_con_yc /= norm_norm;
  norm_con_zc /= norm_norm;
  momentum_factor *= norm_norm;
  double norm_con_x = norm_con_xc + ux * norm_con_tc;
  double norm_con_y = norm_con_yc + uy * norm_con_tc;
  double norm_con_z = norm_con_zc + uz * norm_con_tc;

  // Define unprojected vertical direction in camera frame
  double up_con_xc = 0.0;
  double up_con_yc = 0.0;
  double up_con_zc = 1.0;
  if (image_pole)
  {
    up_con_yc = 1.0;
    up_con_zc = 0.0;
  }

  // Calculate covariant metric in camera frame
  double gcov_xc_xc = gcov(1,1) - u_x / u_t * gcov(1,0) - u_x / u_t * gcov(1,0)
      + u_x * u_x / (u_t * u_t) * gcov(0,0);
  double gcov_xc_yc = gcov(1,2) - u_x / u_t * gcov(2,0) - u_y / u_t * gcov(1,0)
      + u_x * u_y / (u_t * u_t) * gcov(0,0);
  double gcov_xc_zc = gcov(1,3) - u_x / u_t * gcov(3,0) - u_z / u_t * gcov(1,0)
      + u_x * u_z / (u_t * u_t) * gcov(0,0);
  double gcov_yc_yc = gcov(2,2) - u_y / u_t * gcov(2,0) - u_y / u_t * gcov(2,0)
      + u_y * u_y / (u_t * u_t) * gcov(0,0);
  double gcov_yc_zc = gcov(2,3) - u_y / u_t * gcov(3,0) - u_z / u_t * gcov(2,0)
      + u_y * u_z / (u_t * u_t) * gcov(0,0);
  double gcov_zc_zc = gcov(3,3) - u_z / u_t * gcov(3,0) - u_z / u_t * gcov(3,0)
      + u_z * u_z / (u_t * u_t) * gcov(0,0);

  // Calculate camera vertical direction without rotation in camera frame
  double up_norm = up_con_xc * norm_cov_xc + up_con_yc * norm_cov_yc + up_con_zc * norm_cov_zc;
  double vert_con_tc = 0.0;
  double vert_con_xc = up_con_xc - up_norm * norm_con_xc;
  double vert_con_yc = up_con_yc - up_norm * norm_con_yc;
  double vert_con_zc = up_con_zc - up_norm * norm_con_zc;
  double vert_cov_xc =
      gcov_xc_xc * vert_con_xc + gcov_xc_yc * vert_con_yc + gcov_xc_zc * vert_con_zc;
  double vert_cov_yc =
      gcov_xc_yc * vert_con_xc + gcov_yc_yc * vert_con_yc + gcov_yc_zc * vert_con_zc;
  double vert_cov_zc =
      gcov_xc_zc * vert_con_xc + gcov_yc_zc * vert_con_yc + gcov_zc_zc * vert_con_zc;
  double vert_norm =
      std::sqrt(vert_cov_xc * vert_con_xc + vert_cov_yc * vert_con_yc + vert_cov_zc * vert_con_zc);
  vert_cov_xc /= vert_norm;
  vert_cov_yc /= vert_norm;
  vert_cov_zc /= vert_norm;
  vert_con_xc /= vert_norm;
  vert_con_yc /= vert_norm;
  vert_con_zc /= vert_norm;

  // Calculate determinant of metric in camera frame
  double det = gcov_xc_xc * (gcov_yc_yc * gcov_zc_zc - gcov_yc_zc * gcov_yc_zc)
      + gcov_xc_yc * (gcov_yc_zc * gcov_xc_zc - gcov_xc_yc * gcov_zc_zc)
      + gcov_xc_zc * (gcov_xc_yc * gcov_yc_zc - gcov_yc_yc * gcov_xc_zc);
  double det_sqrt = std::sqrt(det);

  // Calculate camera horizontal direction without rotation in camera frame
  double hor_con_tc = 0.0;
  double hor_con_xc = (vert_cov_yc * norm_cov_zc - vert_cov_zc * norm_cov_yc) / det_sqrt;
  double hor_con_yc = (vert_cov_zc * norm_cov_xc - vert_cov_xc * norm_cov_zc) / det_sqrt;
  double hor_con_zc = (vert_cov_xc * norm_cov_yc - vert_cov_yc * norm_cov_xc) / det_sqrt;

  // Calculate camera direction with rotation in camera frame
  double temp_hor_con_xc = hor_con_xc;
  double temp_hor_con_yc = hor_con_yc;
  double temp_hor_con_zc = hor_con_zc;
  double temp_vert_con_xc = vert_con_xc;
  double temp_vert_con_yc = vert_con_yc;
  double temp_vert_con_zc = vert_con_zc;
  hor_con_xc = temp_hor_con_xc * crot - temp_vert_con_xc * srot;
  hor_con_yc = temp_hor_con_yc * crot - temp_vert_con_yc * srot;
  hor_con_zc = temp_hor_con_zc * crot - temp_vert_con_zc * srot;
  vert_con_xc = temp_vert_con_xc * crot + temp_hor_con_xc * srot;
  vert_con_yc = temp_vert_con_yc * crot + temp_hor_con_yc * srot;
  vert_con_zc = temp_vert_con_zc * crot + temp_hor_con_zc * srot;

  // Allocate arrays
  image_position.Allocate(image_resolution, image_resolution, 4);
  image_direction.Allocate(image_resolution, image_resolution, 4);

  // Initialize position and direction based on camera type
  switch (image_camera)
  {
    // Plane with parallel rays
    case Camera::plane:
    {
      #pragma omp parallel for schedule(static)
      for (int m = 0; m < image_resolution; m++)
        for (int l = 0; l < image_resolution; l++)
        {
          // Set pixel position
          double u = (l - image_resolution / 2.0 + 0.5) * bh_m * image_width / image_resolution;
          double v = (m - image_resolution / 2.0 + 0.5) * bh_m * image_width / image_resolution;
          double dtc = u * hor_con_tc + v * vert_con_tc;
          double dxc = u * hor_con_xc + v * vert_con_xc;
          double dyc = u * hor_con_yc + v * vert_con_yc;
          double dzc = u * hor_con_zc + v * vert_con_zc;
          double dt = ut * dtc - (u_x * dxc + u_y * dyc + u_z * dzc) / u_t;
          double dx = dxc + ux * dtc;
          double dy = dyc + uy * dtc;
          double dz = dzc + uz * dtc;
          image_position(m,l,0) = t + dt;
          image_position(m,l,1) = x + dx;
          image_position(m,l,2) = y + dy;
          image_position(m,l,3) = z + dz;

          // Set pixel direction
          image_direction(m,l,1) = norm_con_x;
          image_direction(m,l,2) = norm_con_y;
          image_direction(m,l,3) = norm_con_z;
        }
      break;
    }

    // Point with converging rays
    case Camera::pinhole:
    {
      #pragma omp parallel for schedule(static)
      for (int m = 0; m < image_resolution; m++)
        for (int l = 0; l < image_resolution; l++)
        {
          // Set pixel position
          image_position(m,l,0) = t;
          image_position(m,l,1) = x;
          image_position(m,l,2) = y;
          image_position(m,l,3) = z;

          // Set pixel direction
          double u = (l - image_resolution / 2.0 + 0.5) * bh_m * image_width / image_resolution;
          double v = (m - image_resolution / 2.0 + 0.5) * bh_m * image_width / image_resolution;
          double normalization = std::hypot(u, v, image_r);
          double frac_norm = image_r / normalization;
          double frac_hor = -u / normalization;
          double frac_vert = -v / normalization;
          double dir_con_tc = norm_con_tc;
          double dir_con_xc =
              frac_norm * norm_con_xc + frac_hor * hor_con_xc + frac_vert * vert_con_xc;
          double dir_con_yc =
              frac_norm * norm_con_yc + frac_hor * hor_con_yc + frac_vert * vert_con_yc;
          double dir_con_zc =
              frac_norm * norm_con_zc + frac_hor * hor_con_zc + frac_vert * vert_con_zc;
          double dir_con_x = dir_con_xc + ux * dir_con_tc;
          double dir_con_y = dir_con_yc + uy * dir_con_tc;
          double dir_con_z = dir_con_zc + uz * dir_con_tc;
          image_direction(m,l,1) = dir_con_x;
          image_direction(m,l,2) = dir_con_y;
          image_direction(m,l,3) = dir_con_z;
        }
      break;
    }
  }
  return;
}
