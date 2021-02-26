// Blacklight ray tracer

// C++ headers
#include <algorithm>  // max
#include <cmath>      // abs, acos, cos, sqrt
#include <optional>   // optional

// Blacklight headers
#include "ray_tracer.hpp"
#include "../blacklight.hpp"                   // physics, enums
#include "../athena_reader/athena_reader.hpp"  // AthenaReader
#include "../input_reader/input_reader.hpp"    // InputReader
#include "../utils/array.hpp"                  // Array

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

  // Copy simulation parameters
  if (model_type == ModelType::simulation)
  {
    simulation_coord = p_input_reader->simulation_coord.value();
    simulation_m_msun = p_input_reader->simulation_m_msun.value();
    simulation_rho_cgs = p_input_reader->simulation_rho_cgs.value();
    simulation_interp = p_input_reader->simulation_interp.value();
    if (simulation_interp)
      simulation_block_interp = p_input_reader->simulation_block_interp.value();
  }

  // Copy plasma parameters
  if (model_type == ModelType::simulation)
  {
    plasma_mu = p_input_reader->plasma_mu.value();
    plasma_ne_ni = p_input_reader->plasma_ne_ni.value();
    plasma_model = p_input_reader->plasma_model.value();
    if (plasma_model == PlasmaModel::ti_te_beta)
    {
      plasma_rat_high = p_input_reader->plasma_rat_high.value();
      plasma_rat_low = p_input_reader->plasma_rat_low.value();
    }
    plasma_sigma_max = p_input_reader->plasma_sigma_max.value();
  }

  // Copy fallback parameters
  fallback_nan = p_input_reader->fallback_nan.value();
  if (model_type == ModelType::simulation and not fallback_nan)
  {
    fallback_rho = p_input_reader->fallback_rho.value();
    if (plasma_model == PlasmaModel::ti_te_beta)
      fallback_pgas = p_input_reader->fallback_pgas.value();
    if (plasma_model == PlasmaModel::code_kappa)
      fallback_kappa = p_input_reader->fallback_kappa.value();
  }

  // Copy image parameters
  image_camera = p_input_reader->image_camera.value();
  image_r = p_input_reader->image_r.value();
  image_th = p_input_reader->image_th.value();
  image_ph = p_input_reader->image_ph.value();
  image_urn = p_input_reader->image_urn.value();
  image_uthn = p_input_reader->image_uthn.value();
  image_uphn = p_input_reader->image_uphn.value();
  image_k_r = p_input_reader->image_k_r.value();
  image_k_th = p_input_reader->image_k_th.value();
  image_k_ph = p_input_reader->image_k_ph.value();
  image_rotation = p_input_reader->image_rotation.value();
  image_width = p_input_reader->image_width.value();
  image_resolution = p_input_reader->image_resolution.value();
  image_frequency = p_input_reader->image_frequency.value();
  image_normalization = p_input_reader->image_normalization.value();
  image_pole = p_input_reader->image_pole.value();

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
    if (plasma_model == PlasmaModel::ti_te_beta)
    {
      grid_pgas = p_athena_reader->prim;
      grid_pgas.Slice(5, p_athena_reader->ind_pgas);
    }
    if (plasma_model == PlasmaModel::code_kappa)
    {
      grid_kappa = p_athena_reader->prim;
      grid_kappa.Slice(5, p_athena_reader->ind_kappa);
    }
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
