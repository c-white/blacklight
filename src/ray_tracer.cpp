// Blacklight ray tracer

// C++ headers
#include <algorithm>  // max, min
#include <cmath>      // abs, acos, atan, atan2, cbrt, ceil, cos, cyl_bessel_k, exp, expm1, hypot,
                      // isfinite, pow, sin, sqrt
#include <limits>     // numeric_limits
#include <optional>   // optional
#include <sstream>    // stringstream
#include <string>     // string

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "ray_tracer.hpp"
#include "array.hpp"        // Array
#include "blacklight.hpp"   // math, physics, enumerations
#include "exceptions.hpp"   // BlacklightWarning
#include "read_athena.hpp"  // AthenaReader
#include "read_input.hpp"   // InputReader

//--------------------------------------------------------------------------------------------------

// Ray tracer constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
//   p_athena_reader: pointer to object containing raw simulation data
RayTracer::RayTracer(const InputReader *p_input_reader, const AthenaReader *p_athena_reader)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();

  // Set parameters
  switch (model_type)
  {
    case ModelType::simulation:
    {
      bh_m = 1.0;
      bh_a = p_input_reader->simulation_a.value();
      break;
    }
    case ModelType::formula:
    {
      bh_m = 1.0;
      bh_a = p_input_reader->formula_spin.value();
      break;
    }
  }

  // Copy simulation parameters
  if (model_type == ModelType::simulation)
  {
    simulation_m_msun = p_input_reader->simulation_m_msun.value();
    simulation_rho_cgs = p_input_reader->simulation_rho_cgs.value();
    simulation_coord = p_input_reader->simulation_coord.value();
    simulation_interp = p_input_reader->simulation_interp.value();
    if (simulation_interp)
      simulation_block_interp = p_input_reader->simulation_block_interp.value();
  }

  // Copy plasma parameters
  if (model_type == ModelType::simulation)
  {
    plasma_mu = p_input_reader->plasma_mu.value();
    plasma_ne_ni = p_input_reader->plasma_ne_ni.value();
    plasma_rat_high = p_input_reader->plasma_rat_high.value();
    plasma_rat_low = p_input_reader->plasma_rat_low.value();
    plasma_sigma_max = p_input_reader->plasma_sigma_max.value();
  }

  // Copy formula parameters
  if (model_type == ModelType::formula)
  {
    formula_mass = p_input_reader->formula_mass.value();
    formula_r0 = p_input_reader->formula_r0.value();
    formula_h = p_input_reader->formula_h.value();
    formula_l0 = p_input_reader->formula_l0.value();
    formula_q = p_input_reader->formula_q.value();
    formula_nup = p_input_reader->formula_nup.value();
    formula_cn0 = p_input_reader->formula_cn0.value();
    formula_alpha = p_input_reader->formula_alpha.value();
    formula_a = p_input_reader->formula_a.value();
    formula_beta = p_input_reader->formula_beta.value();
  }

  // Copy fallback parameters
  fallback_nan = p_input_reader->fallback_nan.value();
  if (model_type == ModelType::simulation and not fallback_nan)
  {
    fallback_rho = p_input_reader->fallback_rho.value();
    fallback_pgas = p_input_reader->fallback_pgas.value();
  }

  // Copy image parameters
  im_camera = p_input_reader->im_camera.value();
  im_r = p_input_reader->im_r.value();
  im_th = p_input_reader->im_th.value();
  im_ph = p_input_reader->im_ph.value();
  im_urn = p_input_reader->im_urn.value();
  im_uthn = p_input_reader->im_uthn.value();
  im_uphn = p_input_reader->im_uphn.value();
  im_k_r = p_input_reader->im_k_r.value();
  im_k_th = p_input_reader->im_k_th.value();
  im_k_ph = p_input_reader->im_k_ph.value();
  im_rot = p_input_reader->im_rot.value();
  im_width = p_input_reader->im_width.value();
  im_res = p_input_reader->im_res.value();
  im_freq = p_input_reader->im_freq.value();
  im_norm = p_input_reader->im_norm.value();
  im_pole = p_input_reader->im_pole.value();

  // Copy ray-tracing parameters
  ray_flat = p_input_reader->ray_flat.value();
  ray_terminate = p_input_reader->ray_terminate.value();
  if (ray_terminate != RayTerminate::photon)
    ray_factor = p_input_reader->ray_factor.value();
  ray_step = p_input_reader->ray_step.value();
  ray_max_steps = p_input_reader->ray_max_steps.value();
  ray_max_retries = p_input_reader->ray_max_retries.value();
  ray_tol_abs = p_input_reader->ray_tol_abs.value();
  ray_tol_rel = p_input_reader->ray_tol_rel.value();
  ray_err_factor = p_input_reader->ray_err_factor.value();
  ray_min_factor = p_input_reader->ray_min_factor.value();
  ray_max_factor = p_input_reader->ray_max_factor.value();

  // Copy raw data scalars
  if (model_type == ModelType::simulation and simulation_coord == Coordinates::sph_ks
      and simulation_interp and simulation_block_interp)
    n_3_root = p_athena_reader->n_3_root;

  // Make shallow copies of raw data arrays
  if (model_type == ModelType::simulation)
  {
    if (simulation_interp and simulation_block_interp)
    {
      levels = p_athena_reader->levels;
      locations = p_athena_reader->locations;
    }
    x1f = p_athena_reader->x1f;
    x2f = p_athena_reader->x2f;
    x3f = p_athena_reader->x3f;
    x1v = p_athena_reader->x1v;
    x2v = p_athena_reader->x2v;
    x3v = p_athena_reader->x3v;
    grid_rho = p_athena_reader->prim;
    grid_rho.Slice(5, p_athena_reader->ind_rho);
    grid_pgas = p_athena_reader->prim;
    grid_pgas.Slice(5, p_athena_reader->ind_pgas);
    grid_uu1 = p_athena_reader->prim;
    grid_uu1.Slice(5, p_athena_reader->ind_uu1);
    grid_uu2 = p_athena_reader->prim;
    grid_uu2.Slice(5, p_athena_reader->ind_uu2);
    grid_uu3 = p_athena_reader->prim;
    grid_uu3.Slice(5, p_athena_reader->ind_uu3);
    grid_bb1 = p_athena_reader->bb;
    grid_bb1.Slice(5, p_athena_reader->ind_bb1);
    grid_bb2 = p_athena_reader->bb;
    grid_bb2.Slice(5, p_athena_reader->ind_bb2);
    grid_bb3 = p_athena_reader->bb;
    grid_bb3.Slice(5, p_athena_reader->ind_bb3);
  }

  // Calculate maximum refinement level and number of blocks in x^3-direction at each level
  if (model_type == ModelType::simulation and simulation_coord == Coordinates::sph_ks
      and simulation_interp and simulation_block_interp)
  {
    int n_b = x1f.n2;
    int n_k = x3v.n1;
    max_level = 0;
    for (int b = 0; b < n_b; b++)
      max_level = std::max(max_level, levels(b));
    n_3_level.Allocate(max_level + 1);
    n_3_level(0) = n_3_root / n_k;
    for (int level = 1; level <= max_level; level++)
      n_3_level(level) = n_3_level(level-1) * 2;
  }

  // Calculate black hole mass
  switch (model_type)
  {
    case ModelType::simulation:
    {
      mass_msun = simulation_m_msun;
      break;
    }
    case ModelType::formula:
    {
      mass_msun = formula_mass * physics::c * physics::c / physics::gg_msun;
      break;
    }
  }

  // Calculate termination radii
  double r_horizon = bh_m + std::sqrt(bh_m * bh_m - bh_a * bh_a);
  switch (ray_terminate)
  {
    case RayTerminate::photon:
    {
      r_terminate = 2.0 * bh_m * (1.0 + std::cos(2.0 / 3.0 * std::acos(-std::abs(bh_a) / bh_m)));
      break;
    }
    case RayTerminate::multiplicative:
    {
      r_terminate = r_horizon * ray_factor;
      break;
    }
    case RayTerminate::additive:
    {
      r_terminate = r_horizon + ray_factor;
      break;
    }
  }
}

//--------------------------------------------------------------------------------------------------

// Top-level function for processing raw data into image
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes all data arrays have been set.
void RayTracer::MakeImage()
{
  InitializeCamera();
  InitializeGeodesics();
  IntegrateGeodesics();
  TransformGeodesics();
  switch (model_type)
  {
    case ModelType::simulation:
    {
      SampleSimulationAlongGeodesics();
      IntegrateSimulationRadiation();
      break;
    }
    case ModelType::formula:
    {
      IntegrateFormulaRadiation();
      break;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for setting up camera pixels and initial ray directions
// Inputs: (none)
// Output: (none)
// Notes:
//   Allocates and initializes im_pos and im_dir except for time components of im_dir.
//   Neglects spacetime curvature at camera location.
//   Symbols:
//     n: unit outward normal
//     u: unit right vector
//     v: unit up vector
void RayTracer::InitializeCamera()
{
  // Calculate trigonometric quantities
  double sth = std::sin(im_th);
  double cth = std::cos(im_th);
  double sph = std::sin(im_ph);
  double cph = std::cos(im_ph);
  double srot = std::sin(im_rot);
  double crot = std::cos(im_rot);

  // Calculate camera position
  double t = 0.0;
  double x = sth * (im_r * cph + bh_a * sph);
  double y = sth * (im_r * sph - bh_a * cph);
  double z = im_r * cth;
  if (ray_flat)
  {
    x = im_r * sth * cph;
    y = im_r * sth * sph;
  }

  // Calculate metric in spherical coordinates
  double a2 = bh_a * bh_a;
  double r2 = im_r * im_r;
  double delta = r2 - 2.0 * bh_m * im_r + a2;
  double sigma = r2 + a2 * cth * cth;
  double gcov_r_r = 1.0 + 2.0 * bh_m * im_r / sigma;
  double gcov_r_th = 0.0;
  double gcov_r_ph = -(1.0 + 2.0 * bh_m * im_r / sigma) * bh_a * sth * sth;
  double gcov_th_th = sigma;
  double gcov_th_ph = 0.0;
  double gcov_ph_ph = (r2 + a2 + 2.0 * bh_m * a2 * im_r / sigma * sth * sth) * sth * sth;
  double gcon_t_t = -(1.0 + 2.0 * bh_m * im_r / sigma);
  double gcon_t_r = 2.0 * bh_m * im_r / sigma;
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
  double utn = std::sqrt(1.0 + gcov_r_r * im_urn * im_urn + 2.0 * gcov_r_th * im_urn * im_uthn
      + 2.0 * gcov_r_ph * im_urn * im_uphn + gcov_th_th * im_uthn * im_uthn
      + 2.0 * gcov_th_ph * im_uthn * im_uphn + gcov_ph_ph * im_uphn * im_uphn);
  double ut = utn / alpha;
  double ur = im_urn - beta_con_r / alpha * utn;
  double uth = im_uthn - beta_con_th / alpha * utn;
  double uph = im_uphn - beta_con_ph / alpha * utn;

  // Calculate Jacobian of transformation
  double dx_dr = sth * cph;
  double dy_dr = sth * sph;
  double dz_dr = cth;
  double dx_dth = cth * (im_r * cph + bh_a * sph);
  double dy_dth = cth * (im_r * sph - bh_a * cph);
  double dz_dth = -im_r * sth;
  double dx_dph = sth * (-im_r * sph + bh_a * cph);
  double dy_dph = sth * (im_r * cph + bh_a * sph);
  double dz_dph = 0.0;
  if (ray_flat)
  {
    dx_dth = im_r * cth * cph;
    dy_dth = im_r * cth * sph;
    dx_dph = -im_r * sth * sph;
    dy_dph = im_r * sth * cph;
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
  double k_rn = im_k_r;
  double k_thn = im_k_th;
  double k_phn = im_k_ph;
  double k_tn = -std::sqrt(gcon_rn_rn * k_rn * k_rn + 2.0 * gcon_rn_thn * k_rn * k_thn
      + 2.0 * gcon_rn_phn * k_rn * k_phn + gcon_thn_thn * k_thn * k_thn
      + 2.0 * gcon_thn_phn * k_thn * k_phn + gcon_phn_phn * k_phn * k_phn);
  double k_t = alpha * k_tn + (beta_con_r * k_rn + beta_con_th * k_thn + beta_con_ph * k_phn);

  // Calculate Jacobian of transformation
  double rr2 = x * x + y * y + z * z;
  double dr_dx = im_r * x / (2.0 * r2 - rr2 + a2);
  double dr_dy = im_r * y / (2.0 * r2 - rr2 + a2);
  double dr_dz = (im_r * z + a2 * z / im_r) / (2.0 * r2 - rr2 + a2);
  double dth_dx = z * dr_dx / (r2 * sth);
  double dth_dy = z * dr_dy / (r2 * sth);
  double dth_dz = (z * dr_dz - im_r) / (r2 * sth);
  double dph_dx = -y / (x * x + y * y) - bh_a / (r2 + a2) * dr_dx;
  double dph_dy = x / (x * x + y * y) - bh_a / (r2 + a2) * dr_dy;
  double dph_dz = -bh_a / (r2 + a2) * dr_dz;
  if (ray_flat)
  {
    dr_dx = x / im_r;
    dr_dy = y / im_r;
    dr_dz = z / im_r;
    dth_dx = cth * cph / im_r;
    dth_dy = cth * sph / im_r;
    dth_dz = -sth / im_r;
    dph_dx = -sph / (im_r * sth);
    dph_dy = cph / (im_r * sth);
    dph_dz = 0.0;
  }

  // Calculate photon momentum
  double k_x = dr_dx * im_k_r + dth_dx * im_k_th + dph_dx * im_k_ph;
  double k_y = dr_dy * im_k_r + dth_dy * im_k_th + dph_dy * im_k_ph;
  double k_z = dr_dz * im_k_r + dth_dz * im_k_th + dph_dz * im_k_ph;
  double k_tc = ut * k_t + ux * k_x + uy * k_y + uz * k_z;

  // Calculate momentum normalization
  switch (im_norm)
  {
    case FrequencyNormalization::camera:
    {
      momentum_factor = -im_freq / k_tc;
      break;
    }
    case FrequencyNormalization::infinity:
    {
      momentum_factor = -im_freq / k_t;
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
  if (im_pole)
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
  im_pos.Allocate(im_res, im_res, 4);
  im_dir.Allocate(im_res, im_res, 4);

  // Initialize arrays based on camera type
  switch (im_camera)
  {
    // Plane with parallel rays
    case Camera::plane:
    {
      #pragma omp parallel for schedule(static)
      for (int m = 0; m < im_res; m++)
        for (int l = 0; l < im_res; l++)
        {
          // Set pixel position
          double u = (l - im_res / 2.0 + 0.5) * bh_m * im_width / im_res;
          double v = (m - im_res / 2.0 + 0.5) * bh_m * im_width / im_res;
          double dtc = u * hor_con_tc + v * vert_con_tc;
          double dxc = u * hor_con_xc + v * vert_con_xc;
          double dyc = u * hor_con_yc + v * vert_con_yc;
          double dzc = u * hor_con_zc + v * vert_con_zc;
          double dt = ut * dtc - (u_x * dxc + u_y * dyc + u_z * dzc) / u_t;
          double dx = dxc + ux * dtc;
          double dy = dyc + uy * dtc;
          double dz = dzc + uz * dtc;
          im_pos(m,l,0) = t + dt;
          im_pos(m,l,1) = x + dx;
          im_pos(m,l,2) = y + dy;
          im_pos(m,l,3) = z + dz;

          // Set pixel direction
          im_dir(m,l,1) = norm_con_x;
          im_dir(m,l,2) = norm_con_y;
          im_dir(m,l,3) = norm_con_z;
        }
      break;
    }

    // Point with converging rays
    case Camera::pinhole:
    {
      #pragma omp parallel for schedule(static)
      for (int m = 0; m < im_res; m++)
        for (int l = 0; l < im_res; l++)
        {
          // Set pixel position
          im_pos(m,l,0) = t;
          im_pos(m,l,1) = x;
          im_pos(m,l,2) = y;
          im_pos(m,l,3) = z;

          // Set pixel direction
          double u = (l - im_res / 2.0 + 0.5) * bh_m * im_width / im_res;
          double v = (m - im_res / 2.0 + 0.5) * bh_m * im_width / im_res;
          double normalization = std::hypot(u, v, im_r);
          double frac_norm = im_r / normalization;
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
          im_dir(m,l,1) = dir_con_x;
          im_dir(m,l,2) = dir_con_y;
          im_dir(m,l,3) = dir_con_z;
        }
      break;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for setting up initial conditions for integrating geodesics
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_pos and im_dir have been set except for the time components of im_dir.
//   Initializes time components of im_dir.
//   Lowers all components of im_dir.
//   Quadratic solved as follows:
//     Outside ergosphere: unique positive root.
//     On ergosphere, assuming g_{0i} p^i < 0: unique root, which will be positive
//     Inside ergosphere, assuming g_{0i} p^i < 0: lesser positive root, which remains finite as
//         ergosphere is approached.
void RayTracer::InitializeGeodesics()
{
  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch array
    Array<double> gcov(4, 4);

    // Go through image pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Extract position
        double x[4];
        x[0] = im_pos(m,l,0);
        x[1] = im_pos(m,l,1);
        x[2] = im_pos(m,l,2);
        x[3] = im_pos(m,l,3);

        // Extract spatial components of momentum
        double p[4];
        p[1] = im_dir(m,l,1);
        p[2] = im_dir(m,l,2);
        p[3] = im_dir(m,l,3);

        // Calculate time component of momentum
        CovariantGeodesicMetric(x[1], x[2], x[3], gcov);
        double temp_a = gcov(0,0);
        double temp_b = 0.0;
        for (int a = 1; a < 4; a++)
          temp_b += 2.0 * gcov(0,a) * p[a];
        double temp_c = 0.0;
        for (int a = 1; a < 4; a++)
          for (int b = 1; b < 4; b++)
            temp_c += gcov(a,b) * p[a] * p[b];
        double temp_d = std::sqrt(std::max(temp_b * temp_b - 4.0 * temp_a * temp_c, 0.0));
        p[0] = temp_a == 0 ? -temp_c / (2.0 * temp_b) : (temp_b < 0.0 ?
            2.0 * temp_c / (temp_d - temp_b) : -(temp_b + temp_d) / (2.0 * temp_a));

        // Lower momentum components
        for (int mu = 0; mu < 4; mu++)
        {
          im_dir(m,l,mu) = 0.0;
          for (int nu = 0; nu < 4; nu++)
            im_dir(m,l,mu) += gcov(mu,nu) * p[nu];
        }
      }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating ray positions and directions through space
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_pos and im_dir have been set.
//   Initializes im_steps.
//   Allocates and initializes geodesic_pos, geodesic_dir, geodesic_len, sample_flags, and
//       sample_num.
//   Assumes x^0 is ignorable.
//   Integrates via the Dormand-Prince method (5th-order Runge-Kutta).
//     Method is RK5(4)7M of 1980 JCoAM 6 19.
//     4th-order interpolation follows 1986 MaCom 46 135, which gives a 4th-order estimate of the
//         midpoint and suggests combining this with values and derivatives at both endpoints to fit
//         a unique quartic over the step.
//     See Solving Ordinary Differential Equations I (Hairer, Norsett, Wanner) for coefficients that
//         accomplish the quartic fit without explicitly evaluating the midpoint.
//     All three references have different coefficients for the 4th-order step used in error
//         estimation; the coefficients here follow the original paper.
//     Interpolation is used to take steps small enough to satisfy user input ray_step.
void RayTracer::IntegrateGeodesics()
{
  // Define coefficients
  double a_vals[7][6] = {};
  a_vals[1][0] = 1.0 / 5.0;
  a_vals[2][0] = 3.0 / 40.0;
  a_vals[2][1] = 9.0 / 40.0;
  a_vals[3][0] = 44.0 / 45.0;
  a_vals[3][1] = -56.0 / 15.0;
  a_vals[3][2] = 32.0 / 9.0;
  a_vals[4][0] = 19372.0 / 6561.0;
  a_vals[4][1] = -25360.0 / 2187.0;
  a_vals[4][2] = 64448.0 / 6561.0;
  a_vals[4][3] = -212.0 / 729.0;
  a_vals[5][0] = 9017.0 / 3168.0;
  a_vals[5][1] = -355.0 / 33.0;
  a_vals[5][2] = 46732.0 / 5247.0;
  a_vals[5][3] = 49.0 / 176.0;
  a_vals[5][4] = -5103.0 / 18656.0;
  a_vals[6][0] = 35.0 / 384.0;
  a_vals[6][2] = 500.0 / 1113.0;
  a_vals[6][3] = 125.0 / 192.0;
  a_vals[6][4] = -2187.0 / 6784.0;
  a_vals[6][5] = 11.0 / 84.0;
  double b_vals_5[7] =
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
  double b_vals_4[7] = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0,
      187.0 / 2100.0, 1.0 / 40.0};
  double b_vals_4m[7] = {6025192743.0 / 30085553152.0, 0.0, 51252292925.0 / 65400821598.0,
      -2691868925.0 / 45128329728.0, 187940372067.0 / 1594534317056.0,
      -1776094331.0 / 19743644256.0, 11237099.0 / 235043384.0};
  double d_vals[7] = {-12715105075.0 / 11282082432.0, 0.0, 87487479700.0 / 32700410799.0,
      -10690763975.0 / 1880347072.0, 701980252875.0 / 199316789632.0, -1453857185.0 / 822651844.0,
      69997945.0 / 29380423.0};

  // Define numerical parameter
  double err_power = 0.2;

  // Allocate arrays
  geodesic_pos.Allocate(im_res, im_res, ray_max_steps, 4);
  geodesic_dir.Allocate(im_res, im_res, ray_max_steps, 4);
  geodesic_len.Allocate(im_res, im_res, ray_max_steps);
  geodesic_len.Zero();
  sample_flags.Allocate(im_res, im_res);
  sample_flags.Zero();
  sample_num.Allocate(im_res, im_res);
  sample_num.Zero();

  // Work in parallel
  im_steps = 0;
  int num_bad_geodesics = 0;
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcov(4, 4);
    Array<double> gcon(4, 4);
    Array<double> dgcon(3, 4, 4);
    double y_vals[9];
    double y_vals_temp[9];
    double y_vals_5[9];
    double y_vals_4[9];
    double y_vals_4m[8];
    double k_vals[7][9];
    double r_vals[4][8];

    // Go through image pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Extract initial position
        y_vals[0] = im_pos(m,l,0);
        y_vals[1] = im_pos(m,l,1);
        y_vals[2] = im_pos(m,l,2);
        y_vals[3] = im_pos(m,l,3);

        // Extract initial momentum
        y_vals[4] = im_dir(m,l,0);
        y_vals[5] = im_dir(m,l,1);
        y_vals[6] = im_dir(m,l,2);
        y_vals[7] = im_dir(m,l,3);

        // Set initial proper distance
        y_vals[8] = 0.0;

        // Prepare to take steps
        for (int p = 0; p < 9; p++)
          y_vals_5[p] = y_vals[p];
        double r_new = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);
        double h_new = -r_new * ray_step;
        int num_retry = 0;
        bool previous_fail = false;

        // Take steps
        for (int n = 0; n < ray_max_steps; )
        {
          // Check for too many retries
          if (num_retry > ray_max_retries)
          {
            sample_flags(m,l) = true;
            break;
          }

          // Update step size
          double h = h_new;

          // Copy previous results
          if (not previous_fail and n > 0)
            for (int p = 0; p < 9; p++)
            {
              y_vals[p] = y_vals_5[p];
              k_vals[0][p] = k_vals[6][p];
            }
          if (not previous_fail and n == 0)
            for (int p = 0; p < 9; p++)
              GeodesicSubstep(y_vals, k_vals[0], gcov, gcon, dgcon);
          double r = r_new;
          if (previous_fail)
            r = RadialGeodesicCoordinate(y_vals[1], y_vals[2], y_vals[3]);

          // Calculate substeps
          for (int substep = 1; substep < 7; substep++)
          {
            for (int p = 0; p < 9; p++)
              y_vals_temp[p] = y_vals[p];
            for (int q = 0; q < substep; q++)
              for (int p = 0; p < 9; p++)
                y_vals_temp[p] += a_vals[substep][q] * h * k_vals[q][p];
            GeodesicSubstep(y_vals_temp, k_vals[substep], gcov, gcon, dgcon);
          }

          // Calculate values at end of full step
          for (int p = 0; p < 9; p++)
          {
            y_vals_5[p] = y_vals[p];
            y_vals_4[p] = y_vals[p];
          }
          for (int q = 0; q < 7; q++)
            for (int p = 0; p < 9; p++)
            {
              y_vals_5[p] += b_vals_5[q] * h * k_vals[q][p];
              y_vals_4[p] += b_vals_4[q] * h * k_vals[q][p];
            }
          r_new = RadialGeodesicCoordinate(y_vals_5[1], y_vals_5[2], y_vals_5[3]);

          // Estimate error
          double error = 0.0;
          for (int p = 0; p < 8; p++)
          {
            double y_abs = std::max(std::abs(y_vals[p]), std::abs(y_vals_5[p]));
            double error_scale = ray_tol_abs + ray_tol_rel * y_abs;
            double delta_y = std::abs(y_vals_5[p] - y_vals_4[p]);
            error = std::max(error, delta_y / error_scale);
          }

          // Decide if step is too far
          if (not (error <= 1.0))
          {
            double h_factor = ray_min_factor;
            if (std::isfinite(error))
            {
              double h_factor_ideal = ray_err_factor * std::pow(error, -err_power);
              h_factor = std::max(h_factor_ideal, ray_min_factor);
            }
            h_new = h * h_factor;
            num_retry += 1;
            previous_fail = true;
            continue;
          }
          else
          {
            double h_factor = ray_max_factor;
            if (error > 0.0)
            {
              h_factor = ray_err_factor * std::pow(error, -err_power);
              h_factor = std::max(h_factor, ray_min_factor);
              h_factor = std::min(h_factor, ray_max_factor);
            }
            if (previous_fail)
              h_factor = std::min(h_factor, 1.0);
            h_new = h * h_factor;
            num_retry = 0;
            previous_fail = false;
          }

          // Calculate values at middle of full step
          for (int p = 0; p < 8; p++)
            y_vals_4m[p] = y_vals[p];
          for (int q = 0; q < 7; q++)
            for (int p = 0; p < 8; p++)
              y_vals_4m[p] += b_vals_4m[q] * h * k_vals[q][p];

          // Subdivide full step
          double r_mid = RadialGeodesicCoordinate(y_vals_4m[1], y_vals_4m[2], y_vals_4m[3]);
          double delta_s_step = ray_step * r_mid;
          double delta_s_full = y_vals_5[8] - y_vals[8];
          int num_steps_ideal = static_cast<int>(std::ceil(delta_s_full / delta_s_step));
          delta_s_step = delta_s_full / num_steps_ideal;
          int num_steps_max = ray_max_steps - n;
          int num_steps = num_steps_ideal;
          if (num_steps > num_steps_max)
          {
            num_steps = num_steps_max;
            sample_flags(m,l) = true;
          }

          // Calculate step midpoint if no subdivision necessary
          if (num_steps_ideal == 1)
          {
            geodesic_pos(m,l,n,0) = y_vals_4m[0];
            geodesic_pos(m,l,n,1) = y_vals_4m[1];
            geodesic_pos(m,l,n,2) = y_vals_4m[2];
            geodesic_pos(m,l,n,3) = y_vals_4m[3];
            geodesic_dir(m,l,n,0) = y_vals_4m[4];
            geodesic_dir(m,l,n,1) = y_vals_4m[5];
            geodesic_dir(m,l,n,2) = y_vals_4m[6];
            geodesic_dir(m,l,n,3) = y_vals_4m[7];
            geodesic_len(m,l,n) = h;
          }

          // Calculate interpolating coefficients for subdivisions
          if (num_steps_ideal > 1)
          {
            for (int p = 0; p < 8; p++)
            {
              r_vals[0][p] = y_vals_5[p] - y_vals[p];
              r_vals[1][p] = y_vals[p] - y_vals_5[p] + h * k_vals[0][p];
              r_vals[2][p] = 2.0 * (y_vals_5[p] - y_vals[p]) - h * (k_vals[0][p] + k_vals[6][p]);
              r_vals[3][p] = 0.0;
            }
            for (int q = 0; q < 7; q++)
              for (int p = 0; p < 8; p++)
                r_vals[3][p] += d_vals[q] * h * k_vals[q][p];
          }

          // Calculate subdivided steps
          if (num_steps_ideal > 1)
            for (int nn = 0; nn < num_steps; nn++)
            {
              double frac = (nn + 0.5) / num_steps_ideal;
              for (int p = 0; p < 8; p++)
                y_vals_temp[p] = y_vals[p] + frac * (r_vals[0][p] + (1.0 - frac) * (r_vals[1][p]
                    + frac * (r_vals[2][p] + (1.0 - frac) * r_vals[3][p])));
              geodesic_pos(m,l,n+nn,0) = y_vals_temp[0];
              geodesic_pos(m,l,n+nn,1) = y_vals_temp[1];
              geodesic_pos(m,l,n+nn,2) = y_vals_temp[2];
              geodesic_pos(m,l,n+nn,3) = y_vals_temp[3];
              geodesic_dir(m,l,n+nn,0) = y_vals_temp[4];
              geodesic_dir(m,l,n+nn,1) = y_vals_temp[5];
              geodesic_dir(m,l,n+nn,2) = y_vals_temp[6];
              geodesic_dir(m,l,n+nn,3) = y_vals_temp[7];
              geodesic_len(m,l,n+nn) = h / num_steps_ideal;
            }

          // Renormalize momentum
          ContravariantGeodesicMetric(y_vals_5[1], y_vals_5[2], y_vals_5[3], gcon);
          double temp_a = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              temp_a += gcon(a,b) * y_vals_5[4+a] * y_vals_5[4+b];
          double temp_b = 0.0;
          for (int a = 1; a < 4; a++)
            temp_b += 2.0 * gcon(0,a) * y_vals_5[4] * y_vals_5[4+a];
          double temp_c = gcon(0,0) * y_vals_5[4] * y_vals_5[4];
          double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
          double factor =
              temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
          for (int a = 1; a < 4; a++)
            y_vals_5[4+a] *= factor;

          // Check termination
          sample_num(m,l) += num_steps;
          bool terminate_outer = r_new > im_r and r_new > r;
          bool terminate_inner = r_new < r_terminate;
          if (terminate_outer or terminate_inner)
            break;
          bool last_step = n + num_steps >= ray_max_steps;
          if (last_step)
            sample_flags(m,l) = true;

          // Prepare for next step
          n += num_steps;
        }
      }

    // Renormalize momenta
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        for (int n = 0; n < sample_num(m,l); n++)
        {
          ContravariantGeodesicMetric(geodesic_pos(m,l,n,1), geodesic_pos(m,l,n,2),
              geodesic_pos(m,l,n,3), gcon);
          double temp_a = 0.0;
          for (int a = 1; a < 4; a++)
            for (int b = 1; b < 4; b++)
              temp_a += gcon(a,b) * geodesic_dir(m,l,n,a) * geodesic_dir(m,l,n,b);
          double temp_b = 0.0;
          for (int a = 1; a < 4; a++)
            temp_b += 2.0 * gcon(0,a) * geodesic_dir(m,l,n,0) * geodesic_dir(m,l,n,a);
          double temp_c = gcon(0,0) * geodesic_dir(m,l,n,0) * geodesic_dir(m,l,n,0);
          double temp_d = std::sqrt(temp_b * temp_b - 4.0 * temp_a * temp_c);
          double factor =
              temp_b < 0.0 ? (temp_d - temp_b) / (2.0 * temp_a) : -2.0 * temp_c / (temp_b + temp_d);
          for (int a = 1; a < 4; a++)
            geodesic_dir(m,l,n,a) *= factor;

        }

    // Calculate maximum number of steps actually taken
    #pragma omp for schedule(static) reduction(max:im_steps)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        im_steps = std::max(im_steps, sample_num(m,l));

    // Calculate number of geodesics that do not terminate properly
    #pragma omp for schedule(static) reduction(+:num_bad_geodesics)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        if (sample_flags(m,l))
          num_bad_geodesics++;
  }

  // Report improperly terminated geodesics
  if (num_bad_geodesics > 0)
  {
    std::stringstream message;
    message << num_bad_geodesics << " out of " << im_res * im_res
        << " geodesics terminate unexpectedly.";
    BlacklightWarning(message.str().c_str());
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for transforming geodesics from integrating metric to simulation/formula metric
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_steps, geodesic_pos, geodesic_dir, geodesic_len, and sample_num have been set.
//   Allocates and initializes sample_pos, sample_dir, and sample_len, except transformed from
//       integrating metric to simulation/formula metric and reversed in the sampling dimension.
//   Deallocates geodesic_pos, geodesic_dir, and geodesic_len.
//   Assumes integrating and formula metrics are Cartesian Kerr-Schild.
//   Transformation of time components is trivial.
//   Transformation of geodesic length is trivial.
//   All transformations are trivial if both metrics are Cartesian Kerr-Schild.
void RayTracer::TransformGeodesics()
{
  // Allocate new arrays
  sample_pos.Allocate(im_res, im_res, im_steps, 4);
  sample_dir.Allocate(im_res, im_res, im_steps, 4);
  sample_len.Allocate(im_res, im_res, im_steps);
  sample_len.Zero();

  // Go through samples
  #pragma omp parallel for schedule(static)
  for (int m = 0; m < im_res; m++)
    for (int l = 0; l < im_res; l++)
    {
      int num_steps = sample_num(m,l);
      for (int n = 0; n < num_steps; n++)
      {
        // Skip terminated geodesics
        double len = geodesic_len(m,l,n);
        if (len == 0.0)
          break;

        // Extract Cartesian position
        double t = geodesic_pos(m,l,n,0);
        double x = geodesic_pos(m,l,n,1);
        double y = geodesic_pos(m,l,n,2);
        double z = geodesic_pos(m,l,n,3);

        // Extract Cartesian direction
        double p_t = geodesic_dir(m,l,n,0);
        double p_x = geodesic_dir(m,l,n,1);
        double p_y = geodesic_dir(m,l,n,2);
        double p_z = geodesic_dir(m,l,n,3);

        // Prepare to find positions and momenta
        double x1, x2, x3;
        double p_1, p_2, p_3;

        // Account for model
        switch (model_type)
        {
          // Simulation output
          case ModelType::simulation:
          default:
          {
            // Account for simulation metric
            switch (simulation_coord)
            {
              // Spherical Kerr-Schild
              case Coordinates::sph_ks:
              default:
              {
                // Calculate spherical position
                double a2 = bh_a * bh_a;
                double rr2 = x * x + y * y + z * z;
                double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
                double r = std::sqrt(r2);
                double th = std::acos(z / r);
                double ph = std::atan2(y, x) + std::atan(bh_a / r);
                ph += ph < 0.0 ? 2.0 * math::pi : 0.0;
                ph -= ph > 2.0 * math::pi ? 2.0 * math::pi : 0.0;
                double sth = std::sin(th);
                double cth = std::cos(th);
                double sph = std::sin(ph);
                double cph = std::cos(ph);

                // Calculate Jacobian of transformation
                double dx_dr = sth * cph;
                double dy_dr = sth * sph;
                double dz_dr = cth;
                double dx_dth = cth * (r * cph + bh_a * sph);
                double dy_dth = cth * (r * sph - bh_a * cph);
                double dz_dth = -r * sth;
                double dx_dph = sth * (-r * sph + bh_a * cph);
                double dy_dph = sth * (r * cph + bh_a * sph);
                double dz_dph = 0.0;

                // Calculate spherical direction
                double p_r = dx_dr * p_x + dy_dr * p_y + dz_dr * p_z;
                double p_th = dx_dth * p_x + dy_dth * p_y + dz_dth * p_z;
                double p_ph = dx_dph * p_x + dy_dph * p_y + dz_dph * p_z;

                // Assign position and direction
                x1 = r;
                x2 = th;
                x3 = ph;
                p_1 = p_r;
                p_2 = p_th;
                p_3 = p_ph;
                break;
              }

              // Cartesian Kerr-Schild
              case Coordinates::cart_ks:
              {
                x1 = x;
                x2 = y;
                x3 = z;
                p_1 = p_x;
                p_2 = p_y;
                p_3 = p_z;
                break;
              }
            }
            break;
          }

          // Formula
          case ModelType::formula:
          {
            x1 = x;
            x2 = y;
            x3 = z;
            p_1 = p_x;
            p_2 = p_y;
            p_3 = p_z;
            break;
          }
        }

        // Set new arrays in reverse order
        sample_pos(m,l,num_steps-1-n,0) = t;
        sample_pos(m,l,num_steps-1-n,1) = x1;
        sample_pos(m,l,num_steps-1-n,2) = x2;
        sample_pos(m,l,num_steps-1-n,3) = x3;
        sample_dir(m,l,num_steps-1-n,0) = p_t;
        sample_dir(m,l,num_steps-1-n,1) = p_1;
        sample_dir(m,l,num_steps-1-n,2) = p_2;
        sample_dir(m,l,num_steps-1-n,3) = p_3;
        sample_len(m,l,num_steps-1-n) = -len;
      }
    }

  // Deallocate old arrays
  geodesic_pos.Deallocate();
  geodesic_dir.Deallocate();
  geodesic_len.Deallocate();
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for resampling simulation cell data onto rays.
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes im_steps, sample_flags, sample_num, sample_pos, and sample_len have been set.
//   Allocates and initializes sample_rho, sample_pgas, sample_uu1, sample_uu2, sample_uu3,
//       sample_bb1, sample_bb2, and sample_bb3.
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
  sample_rho.Allocate(im_res, im_res, im_steps);
  sample_pgas.Allocate(im_res, im_res, im_steps);
  sample_uu1.Allocate(im_res, im_res, im_steps);
  sample_uu2.Allocate(im_res, im_res, im_steps);
  sample_uu3.Allocate(im_res, im_res, im_steps);
  sample_bb1.Allocate(im_res, im_res, im_steps);
  sample_bb2.Allocate(im_res, im_res, im_steps);
  sample_bb3.Allocate(im_res, im_res, im_steps);

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
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Extract number of steps along this geodesic
        int num_steps = sample_num(m,l);

        // Set NaN fallback values if geodesic poorly terminated
        if (fallback_nan and sample_flags(m,l))
        {
          for (int n = 0; n < num_steps; n++)
          {
            sample_rho(m,l,n) = std::numeric_limits<float>::quiet_NaN();
            sample_pgas(m,l,n) = std::numeric_limits<float>::quiet_NaN();
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
                sample_pgas(m,l,n) = std::numeric_limits<float>::quiet_NaN();
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
                sample_pgas(m,l,n) = fallback_pgas;
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
            sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
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
            if (sample_pgas(m,l,n) <= 0.0f)
              sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
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
            sample_pgas(m,l,n) =
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
            if (sample_pgas(m,l,n) <= 0.0f)
              sample_pgas(m,l,n) = grid_pgas(b,k,j,i);
          }
        }
      }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on simulation data
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_dir, sample_len, sample_rho, sample_pgas, sample_uu1, sample_uu2,
//       sample_uu3, sample_bb1, sample_bb2, and sample_bb3 have been set.
//   Allocates and initializes image.
//   Assumes x^0 is ignorable.
//   References symphony paper 2016 ApJ 822 34 (S).
void RayTracer::IntegrateSimulationRadiation()
{
  // Allocate image array
  image.Allocate(im_res, im_res);
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
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
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
          double pgas = sample_pgas(m,l,n);
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
          double beta_inv = b_sq / (2.0 * pgas);
          double tt_rat = (plasma_rat_high + plasma_rat_low * beta_inv * beta_inv)
              / (1.0 + beta_inv * beta_inv);
          double kb_tt_tot_cgs = plasma_mu * physics::m_p * physics::c * physics::c * pgas / rho;
          double kb_tt_e_cgs = (plasma_ne_ni + 1.0) / (plasma_ne_ni + tt_rat) * kb_tt_tot_cgs;
          double theta_e = kb_tt_e_cgs / (physics::m_e * physics::c * physics::c);
          double nu_s_cgs = 2.0 / 9.0 * nu_c_cgs * theta_e * theta_e * sin_theta;

          // Calculate emission coefficient in CGS units (S 24), skipping contribution if necessary
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
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        image(m,l) *= im_freq * im_freq * im_freq;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for integrating radiative transfer equation based on formula
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes sample_num, sample_dir, and sample_len have been set.
//   Allocates and initializes image.
//   Assumes x^0 is ignorable.
//   References code comparison paper 2020 ApJ 897 148 (C).
void RayTracer::IntegrateFormulaRadiation()
{
  // Allocate image array
  image.Allocate(im_res, im_res);
  image.Zero();

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch arrays
    Array<double> gcon(4, 4);

    // Go through pixels
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
      {
        // Set pixel to NaN if ray has problem
        if (fallback_nan and sample_flags(m,l))
        {
          image(m,l) = std::numeric_limits<double>::quiet_NaN();
          continue;
        }

        // Go through samples
        int num_steps = sample_num(m,l);
        for (int n = 0; n < num_steps; n++)
        {
          // Extract geodesic position and momentum
          double x = sample_pos(m,l,n,1);
          double y = sample_pos(m,l,n,2);
          double z = sample_pos(m,l,n,3);
          double p_0 = sample_dir(m,l,n,0);
          double p_1 = sample_dir(m,l,n,1);
          double p_2 = sample_dir(m,l,n,2);
          double p_3 = sample_dir(m,l,n,3);

          // Calculate coordinates
          double r = RadialGeodesicCoordinate(x, y, z);
          double rr = std::sqrt(r * r - z * z);
          double cth = z / r;
          double sth = std::sqrt(1.0 - cth * cth);
          double ph = std::atan2(y, x) + std::atan(bh_a / r);
          double sph = std::sin(ph);
          double cph = std::cos(ph);

          // Calculate metric
          ContravariantGeodesicMetric(x, y, z, gcon);
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
          double u1 = sth * cph * ur + cth * (r * cph + bh_a * sph) * uth
              + sth * (-r * sph + bh_a * cph) * uph;
          double u2 = sth * sph * ur + cth * (r * sph - bh_a * cph) * uth
              + sth * (r * cph + bh_a * sph) * uph;
          double u3 = cth * ur - r * sth * uth;

          // Calculate fluid-frame number density (C 5)
          double n_n0_fluid = std::exp(-0.5
              * (r * r / (formula_r0 * formula_r0) + formula_h * formula_h * cth * cth));

          // Calculate frequency in CGS units
          double nu_fluid_cgs = -(u0 * p_0 + u1 * p_1 + u2 * p_2 + u3 * p_3) * momentum_factor;

          // Calculate emission coefficient in CGS units (C 9-10)
          double j_nu_fluid_cgs =
              formula_cn0 * n_n0_fluid * std::pow(nu_fluid_cgs / formula_nup, -formula_alpha);
          double j_nu_invariant_cgs = j_nu_fluid_cgs / (nu_fluid_cgs * nu_fluid_cgs);

          // Calculate absorption coefficient in CGS units (C 11-12)
          double k_nu_fluid_cgs = formula_a * formula_cn0 * n_n0_fluid
              * std::pow(nu_fluid_cgs / formula_nup, -formula_beta - formula_alpha);
          double k_nu_invariant_cgs = k_nu_fluid_cgs * nu_fluid_cgs;

          // Calculate change in invariant intensity
          double delta_lambda = sample_len(m,l,n);
          double delta_lambda_cgs = delta_lambda * formula_mass / momentum_factor;
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
    #pragma omp for schedule(static)
    for (int m = 0; m < im_res; m++)
      for (int l = 0; l < im_res; l++)
        image(m,l) *= im_freq * im_freq * im_freq;
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating radial coordinate given location in coordinates used for geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   returned value: radial coordinate
// Notes:
//   Assumes Cartesian Kerr-Schild coordinates.
double RayTracer::RadialGeodesicCoordinate(double x, double y, double z)
{
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  return r;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating covariant metric components for integrating geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   gcov: components set
// Notes:
//   Assumes gcov is allocated to be 4*4.
//   Assumes Cartesian Kerr-Schild coordinates (assumes Minkowski coordinates if ray_flat == true).
void RayTracer::CovariantGeodesicMetric(double x, double y, double z, Array<double> &gcov)
{
  // Handle flat case
  if (ray_flat)
  {
    gcov.Zero();
    gcov(0,0) = -1.0;
    gcov(1,1) = 1.0;
    gcov(2,2) = 1.0;
    gcov(3,3) = 1.0;
    return;
  }

  // Calculate useful quantities
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

  // Calculate null vector
  double l_0 = 1.0;
  double l_1 = (r * x + bh_a * y) / (r2 + a2);
  double l_2 = (r * y - bh_a * x) / (r2 + a2);
  double l_3 = z / r;

  // Calculate metric components
  gcov(0,0) = f * l_0 * l_0 - 1.0;
  gcov(0,1) = f * l_0 * l_1;
  gcov(0,2) = f * l_0 * l_2;
  gcov(0,3) = f * l_0 * l_3;
  gcov(1,0) = f * l_1 * l_0;
  gcov(1,1) = f * l_1 * l_1 + 1.0;
  gcov(1,2) = f * l_1 * l_2;
  gcov(1,3) = f * l_1 * l_3;
  gcov(2,0) = f * l_2 * l_0;
  gcov(2,1) = f * l_2 * l_1;
  gcov(2,2) = f * l_2 * l_2 + 1.0;
  gcov(2,3) = f * l_2 * l_3;
  gcov(3,0) = f * l_3 * l_0;
  gcov(3,1) = f * l_3 * l_1;
  gcov(3,2) = f * l_3 * l_2;
  gcov(3,3) = f * l_3 * l_3 + 1.0;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating contravariant metric components for integrating geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   gcon: components set
// Notes:
//   Assumes gcon is allocated to be 4*4.
//   Assumes Cartesian Kerr-Schild coordinates (assumes Minkowski coordinates if ray_flat == true).
void RayTracer::ContravariantGeodesicMetric(double x, double y, double z, Array<double> &gcon)
{
  // Handle flat case
  if (ray_flat)
  {
    gcon.Zero();
    gcon(0,0) = -1.0;
    gcon(1,1) = 1.0;
    gcon(2,2) = 1.0;
    gcon(3,3) = 1.0;
    return;
  }

  // Calculate useful quantities
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

  // Calculate null vector
  double l0 = -1.0;
  double l1 = (r * x + bh_a * y) / (r2 + a2);
  double l2 = (r * y - bh_a * x) / (r2 + a2);
  double l3 = z / r;

  // Calculate metric components
  gcon(0,0) = -f * l0 * l0 - 1.0;
  gcon(0,1) = -f * l0 * l1;
  gcon(0,2) = -f * l0 * l2;
  gcon(0,3) = -f * l0 * l3;
  gcon(1,0) = -f * l1 * l0;
  gcon(1,1) = -f * l1 * l1 + 1.0;
  gcon(1,2) = -f * l1 * l2;
  gcon(1,3) = -f * l1 * l3;
  gcon(2,0) = -f * l2 * l0;
  gcon(2,1) = -f * l2 * l1;
  gcon(2,2) = -f * l2 * l2 + 1.0;
  gcon(2,3) = -f * l2 * l3;
  gcon(3,0) = -f * l3 * l0;
  gcon(3,1) = -f * l3 * l1;
  gcon(3,2) = -f * l3 * l2;
  gcon(3,3) = -f * l3 * l3 + 1.0;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating derivatives of contravariant metric components for integrating geodesics
// Inputs:
//   x, y, z: coordinates
// Output:
//   dgcon: components set
// Notes:
//   Assumes dgcon is allocated to be 3*4*4.
//   Assumes Cartesian Kerr-Schild coordinates (assumes Minkowski coordinates if ray_flat == true).
void RayTracer::ContravariantGeodesicMetricDerivative(double x, double y, double z,
    Array<double> &dgcon)
{
  // Handle flat case
  if (ray_flat)
  {
    dgcon.Zero();
    return;
  }

  // Calculate useful quantities
  double a2 = bh_a * bh_a;
  double rr2 = x * x + y * y + z * z;
  double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
  double r = std::sqrt(r2);
  double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

  // Calculate null vector
  double l0 = -1.0;
  double l1 = (r * x + bh_a * y) / (r2 + a2);
  double l2 = (r * y - bh_a * x) / (r2 + a2);
  double l3 = z / r;

  // Calculate scalar derivatives
  double dr_dx = r * x / (2.0 * r2 - rr2 + a2);
  double dr_dy = r * y / (2.0 * r2 - rr2 + a2);
  double dr_dz = (r * z + a2 * z / r) / (2.0 * r2 - rr2 + a2);
  double df_dx = -(r2 * r2 - 3.0 * a2 * z * z) * dr_dx / (r * (r2 * r2 + a2 * z * z)) * f;
  double df_dy = -(r2 * r2 - 3.0 * a2 * z * z) * dr_dy / (r * (r2 * r2 + a2 * z * z)) * f;
  double df_dz =
      -((r2 * r2 - 3.0 * a2 * z * z) * dr_dz + 2.0 * a2 * r * z) / (r * (r2 * r2 + a2 * z * z)) * f;

  // Calculate vector derivatives
  double dl0_dx = 0.0;
  double dl0_dy = 0.0;
  double dl0_dz = 0.0;
  double dl1_dx = ((x - 2.0 * r * l1) * dr_dx + r) / (r2 + a2);
  double dl1_dy = ((x - 2.0 * r * l1) * dr_dy + bh_a) / (r2 + a2);
  double dl1_dz = (x - 2.0 * r * l1) * dr_dz / (r2 + a2);
  double dl2_dx = ((y - 2.0 * r * l2) * dr_dx - bh_a) / (r2 + a2);
  double dl2_dy = ((y - 2.0 * r * l2) * dr_dy + r) / (r2 + a2);
  double dl2_dz = (y - 2.0 * r * l2) * dr_dz / (r2 + a2);
  double dl3_dx = -z / r2 * dr_dx;
  double dl3_dy = -z / r2 * dr_dy;
  double dl3_dz = -z / r2 * dr_dz + 1.0 / r;

  // Calculate metric component x-derivatives
  dgcon(0,0,0) = -(df_dx * l0 * l0 + f * dl0_dx * l0 + f * l0 * dl0_dx);
  dgcon(0,0,1) = -(df_dx * l0 * l1 + f * dl0_dx * l1 + f * l0 * dl1_dx);
  dgcon(0,0,2) = -(df_dx * l0 * l2 + f * dl0_dx * l2 + f * l0 * dl2_dx);
  dgcon(0,0,3) = -(df_dx * l0 * l3 + f * dl0_dx * l3 + f * l0 * dl3_dx);
  dgcon(0,1,0) = -(df_dx * l1 * l0 + f * dl1_dx * l0 + f * l1 * dl0_dx);
  dgcon(0,1,1) = -(df_dx * l1 * l1 + f * dl1_dx * l1 + f * l1 * dl1_dx);
  dgcon(0,1,2) = -(df_dx * l1 * l2 + f * dl1_dx * l2 + f * l1 * dl2_dx);
  dgcon(0,1,3) = -(df_dx * l1 * l3 + f * dl1_dx * l3 + f * l1 * dl3_dx);
  dgcon(0,2,0) = -(df_dx * l2 * l0 + f * dl2_dx * l0 + f * l2 * dl0_dx);
  dgcon(0,2,1) = -(df_dx * l2 * l1 + f * dl2_dx * l1 + f * l2 * dl1_dx);
  dgcon(0,2,2) = -(df_dx * l2 * l2 + f * dl2_dx * l2 + f * l2 * dl2_dx);
  dgcon(0,2,3) = -(df_dx * l2 * l3 + f * dl2_dx * l3 + f * l2 * dl3_dx);
  dgcon(0,3,0) = -(df_dx * l3 * l0 + f * dl3_dx * l0 + f * l3 * dl0_dx);
  dgcon(0,3,1) = -(df_dx * l3 * l1 + f * dl3_dx * l1 + f * l3 * dl1_dx);
  dgcon(0,3,2) = -(df_dx * l3 * l2 + f * dl3_dx * l2 + f * l3 * dl2_dx);
  dgcon(0,3,3) = -(df_dx * l3 * l3 + f * dl3_dx * l3 + f * l3 * dl3_dx);

  // Calculate metric component y-derivatives
  dgcon(1,0,0) = -(df_dy * l0 * l0 + f * dl0_dy * l0 + f * l0 * dl0_dy);
  dgcon(1,0,1) = -(df_dy * l0 * l1 + f * dl0_dy * l1 + f * l0 * dl1_dy);
  dgcon(1,0,2) = -(df_dy * l0 * l2 + f * dl0_dy * l2 + f * l0 * dl2_dy);
  dgcon(1,0,3) = -(df_dy * l0 * l3 + f * dl0_dy * l3 + f * l0 * dl3_dy);
  dgcon(1,1,0) = -(df_dy * l1 * l0 + f * dl1_dy * l0 + f * l1 * dl0_dy);
  dgcon(1,1,1) = -(df_dy * l1 * l1 + f * dl1_dy * l1 + f * l1 * dl1_dy);
  dgcon(1,1,2) = -(df_dy * l1 * l2 + f * dl1_dy * l2 + f * l1 * dl2_dy);
  dgcon(1,1,3) = -(df_dy * l1 * l3 + f * dl1_dy * l3 + f * l1 * dl3_dy);
  dgcon(1,2,0) = -(df_dy * l2 * l0 + f * dl2_dy * l0 + f * l2 * dl0_dy);
  dgcon(1,2,1) = -(df_dy * l2 * l1 + f * dl2_dy * l1 + f * l2 * dl1_dy);
  dgcon(1,2,2) = -(df_dy * l2 * l2 + f * dl2_dy * l2 + f * l2 * dl2_dy);
  dgcon(1,2,3) = -(df_dy * l2 * l3 + f * dl2_dy * l3 + f * l2 * dl3_dy);
  dgcon(1,3,0) = -(df_dy * l3 * l0 + f * dl3_dy * l0 + f * l3 * dl0_dy);
  dgcon(1,3,1) = -(df_dy * l3 * l1 + f * dl3_dy * l1 + f * l3 * dl1_dy);
  dgcon(1,3,2) = -(df_dy * l3 * l2 + f * dl3_dy * l2 + f * l3 * dl2_dy);
  dgcon(1,3,3) = -(df_dy * l3 * l3 + f * dl3_dy * l3 + f * l3 * dl3_dy);

  // Calculate metric component z-derivatives
  dgcon(2,0,0) = -(df_dz * l0 * l0 + f * dl0_dz * l0 + f * l0 * dl0_dz);
  dgcon(2,0,1) = -(df_dz * l0 * l1 + f * dl0_dz * l1 + f * l0 * dl1_dz);
  dgcon(2,0,2) = -(df_dz * l0 * l2 + f * dl0_dz * l2 + f * l0 * dl2_dz);
  dgcon(2,0,3) = -(df_dz * l0 * l3 + f * dl0_dz * l3 + f * l0 * dl3_dz);
  dgcon(2,1,0) = -(df_dz * l1 * l0 + f * dl1_dz * l0 + f * l1 * dl0_dz);
  dgcon(2,1,1) = -(df_dz * l1 * l1 + f * dl1_dz * l1 + f * l1 * dl1_dz);
  dgcon(2,1,2) = -(df_dz * l1 * l2 + f * dl1_dz * l2 + f * l1 * dl2_dz);
  dgcon(2,1,3) = -(df_dz * l1 * l3 + f * dl1_dz * l3 + f * l1 * dl3_dz);
  dgcon(2,2,0) = -(df_dz * l2 * l0 + f * dl2_dz * l0 + f * l2 * dl0_dz);
  dgcon(2,2,1) = -(df_dz * l2 * l1 + f * dl2_dz * l1 + f * l2 * dl1_dz);
  dgcon(2,2,2) = -(df_dz * l2 * l2 + f * dl2_dz * l2 + f * l2 * dl2_dz);
  dgcon(2,2,3) = -(df_dz * l2 * l3 + f * dl2_dz * l3 + f * l2 * dl3_dz);
  dgcon(2,3,0) = -(df_dz * l3 * l0 + f * dl3_dz * l0 + f * l3 * dl0_dz);
  dgcon(2,3,1) = -(df_dz * l3 * l1 + f * dl3_dz * l1 + f * l3 * dl1_dz);
  dgcon(2,3,2) = -(df_dz * l3 * l2 + f * dl3_dz * l2 + f * l3 * dl2_dz);
  dgcon(2,3,3) = -(df_dz * l3 * l3 + f * dl3_dz * l3 + f * l3 * dl3_dz);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating covariant metric components in simulation coordinates
// Inputs:
//   x1, x2, x3: coordinates
// Output:
//   gcov: components set
// Notes:
//   Assumes gcov is allocated to be 4*4.
void RayTracer::CovariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcov)
{
  // Account for simulation metric
  switch (simulation_coord)
  {
    // Spherical Kerr-Schild
    case Coordinates::sph_ks:
    {
      // Calculate useful quantities
      double r = x1;
      double th = x2;
      double sth = std::sin(th);
      double cth = std::cos(th);
      double sigma = r * r + bh_a * bh_a * cth * cth;

      // Calculate metric components
      gcov(0,0) = -(1.0 - 2.0 * bh_m * r / sigma);
      gcov(0,1) = 2.0 * bh_m * r / sigma;
      gcov(0,2) = 0.0;
      gcov(0,3) = -2.0 * bh_m * bh_a * r * sth * sth / sigma;
      gcov(1,0) = 2.0 * bh_m * r / sigma;
      gcov(1,1) = 1.0 + 2.0 * bh_m * r / sigma;
      gcov(1,2) = 0.0;
      gcov(1,3) = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * sth * sth;
      gcov(2,0) = 0.0;
      gcov(2,1) = 0.0;
      gcov(2,2) = sigma;
      gcov(2,3) = 0.0;
      gcov(3,0) = -2.0 * bh_m * bh_a * r * sth * sth / sigma;
      gcov(3,1) = -(1.0 + 2.0 * bh_m * r / sigma) * bh_a * sth * sth;
      gcov(3,2) = 0.0;
      gcov(3,3) =
          (r * r + bh_a * bh_a + 2.0 * bh_m * bh_a * bh_a * r * sth * sth / sigma) * sth * sth;
      break;
    }

    // Cartesian Kerr-Schild
    case Coordinates::cart_ks:
    {
      // Calculate useful quantities
      double x = x1;
      double y = x2;
      double z = x3;
      double a2 = bh_a * bh_a;
      double rr2 = x * x + y * y + z * z;
      double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
      double r = std::sqrt(r2);
      double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

      // Calculate null vector
      double l_0 = 1.0;
      double l_1 = (r * x + bh_a * y) / (r2 + a2);
      double l_2 = (r * y - bh_a * x) / (r2 + a2);
      double l_3 = z / r;

      // Calculate metric components
      gcov(0,0) = f * l_0 * l_0 - 1.0;
      gcov(0,1) = f * l_0 * l_1;
      gcov(0,2) = f * l_0 * l_2;
      gcov(0,3) = f * l_0 * l_3;
      gcov(1,0) = f * l_1 * l_0;
      gcov(1,1) = f * l_1 * l_1 + 1.0;
      gcov(1,2) = f * l_1 * l_2;
      gcov(1,3) = f * l_1 * l_3;
      gcov(2,0) = f * l_2 * l_0;
      gcov(2,1) = f * l_2 * l_1;
      gcov(2,2) = f * l_2 * l_2 + 1.0;
      gcov(2,3) = f * l_2 * l_3;
      gcov(3,0) = f * l_3 * l_0;
      gcov(3,1) = f * l_3 * l_1;
      gcov(3,2) = f * l_3 * l_2;
      gcov(3,3) = f * l_3 * l_3 + 1.0;
      break;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for calculating contravariant metric components in simulation coordinates
// Inputs:
//   x1, x2, x3: coordinates
// Output:
//   gcon: components set
// Notes:
//   Assumes gcon is allocated to be 4*4.
void RayTracer::ContravariantCoordinateMetric(double x1, double x2, double x3, Array<double> &gcon)
{
  // Account for simulation metric
  switch (simulation_coord)
  {
    // Spherical Kerr-Schild
    case Coordinates::sph_ks:
    {
      // Calculate useful quantities
      double r = x1;
      double th = x2;
      double sth = std::sin(th);
      double cth = std::cos(th);
      double delta = r * r - 2.0 * bh_m * r + bh_a * bh_a;
      double sigma = r * r + bh_a * bh_a * cth * cth;

      // Calculate metric components
      gcon(0,0) = -(1.0 + 2.0 * bh_m * r / sigma);
      gcon(0,1) = 2.0 * bh_m * r / sigma;
      gcon(0,2) = 0.0;
      gcon(0,3) = 0.0;
      gcon(1,0) = 2.0 * bh_m * r / sigma;
      gcon(1,1) = delta / sigma;
      gcon(1,2) = 0.0;
      gcon(1,3) = bh_a / sigma;
      gcon(2,0) = 0.0;
      gcon(2,1) = 0.0;
      gcon(2,2) = 1.0 / sigma;
      gcon(2,3) = 0.0;
      gcon(3,0) = 0.0;
      gcon(3,1) = bh_a / sigma;
      gcon(3,2) = 0.0;
      gcon(3,3) = 1.0 / (sigma * sth * sth);
      break;
    }

    // Cartesian Kerr-Schild
    case Coordinates::cart_ks:
    {
      // Calculate useful quantities
      double x = x1;
      double y = x2;
      double z = x3;
      double a2 = bh_a * bh_a;
      double rr2 = x * x + y * y + z * z;
      double r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * bh_a * z));
      double r = std::sqrt(r2);
      double f = 2.0 * bh_m * r * r * r / (r2 * r2 + a2 * z * z);

      // Calculate null vector
      double l0 = -1.0;
      double l1 = (r * x + bh_a * y) / (r2 + a2);
      double l2 = (r * y - bh_a * x) / (r2 + a2);
      double l3 = z / r;

      // Calculate metric components
      gcon(0,0) = -f * l0 * l0 - 1.0;
      gcon(0,1) = -f * l0 * l1;
      gcon(0,2) = -f * l0 * l2;
      gcon(0,3) = -f * l0 * l3;
      gcon(1,0) = -f * l1 * l0;
      gcon(1,1) = -f * l1 * l1 + 1.0;
      gcon(1,2) = -f * l1 * l2;
      gcon(1,3) = -f * l1 * l3;
      gcon(2,0) = -f * l2 * l0;
      gcon(2,1) = -f * l2 * l1;
      gcon(2,2) = -f * l2 * l2 + 1.0;
      gcon(2,3) = -f * l2 * l3;
      gcon(3,0) = -f * l3 * l0;
      gcon(3,1) = -f * l3 * l1;
      gcon(3,2) = -f * l3 * l2;
      gcon(3,3) = -f * l3 * l3 + 1.0;
      break;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for taking single forward-Euler substep in time
// Inputs:
//   y: dependent variables (positions, momenta, proper distace)
// Outputs:
//   k: derivatives with respect to independent variable (affine parameter)
//   gcov: components set
//   gcon: components set
//   dgcon: components set
// Notes:
//   Integrates following equations:
//     d(x^mu) / d(lambda) = g^{mu nu} p_nu
//     d(p_0) / d(lambda) = 0
//     d(p_i) / d(lambda) = -1/2 * d(g^{mu nu}) / d(x^i) p_mu p_nu
//     d(s) / d(lambda) = (g_{i j} g^{i mu} g^{j nu} p_mu p_nu)^(1/2)
//   Assumes x^0 is ignorable.
void RayTracer::GeodesicSubstep(double y[9], double k[9], Array<double> &gcov, Array<double> &gcon,
    Array<double> &dgcon)
{
  CovariantGeodesicMetric(y[1], y[2], y[3], gcov);
  ContravariantGeodesicMetric(y[1], y[2], y[3], gcon);
  ContravariantGeodesicMetricDerivative(y[1], y[2], y[3], dgcon);
  for (int p = 0; p < 9; p++)
    k[p] = 0.0;
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      k[mu] += gcon(mu,nu) * y[4+nu];
  for (int a = 1; a < 4; a++)
    for (int mu = 0; mu < 4; mu++)
      for (int nu = 0; nu < 4; nu++)
        k[4+a] -= 0.5 * dgcon(a-1,mu,nu) * y[4+mu] * y[4+nu];
  for (int a = 1; a < 4; a++)
    for (int b = 1; b < 4; b++)
      for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
          k[8] += gcov(a,b) * gcon(a,mu) * gcon(b,nu) * y[4+mu] * y[4+nu];
  k[8] = -std::sqrt(k[8]);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for finding cell values in a given or nearby block
// Inputs:
//   b: block index
//   k, j, i: cell indices for x^3, x^2, and x^1, possibly 1 beyond valid range
// Outputs:
//   vals: density, gas pressure, velocity, and magnetic field set
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
    vals[1] = grid_pgas(b,k,j,i);
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
      vals[1] = grid_pgas(b_alt,k_sought,j_sought,i_sought);
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
        vals[1] = grid_pgas(b_alt,k_sought,j_sought,i_sought);
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
      vals[1] = (static_cast<double>(grid_pgas(b_alt,k_sought,j_sought,i_sought))
          + static_cast<double>(grid_pgas(b_alt,k_sought,j_sought,i_sought+1))
          + static_cast<double>(grid_pgas(b_alt,k_sought,j_sought+1,i_sought))
          + static_cast<double>(grid_pgas(b_alt,k_sought,j_sought+1,i_sought+1))
          + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought,i_sought))
          + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought,i_sought+1))
          + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought+1,i_sought))
          + static_cast<double>(grid_pgas(b_alt,k_sought+1,j_sought+1,i_sought+1))) / 8.0;
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
