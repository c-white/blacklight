// Blacklight radiation integrator

// C++ headers
#include <algorithm>  // max
#include <optional>   // optional

// Library headers
#include <omp.h>  // omp_get_wtime

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"                               // physics, enums
#include "../athena_reader/athena_reader.hpp"              // AthenaReader
#include "../geodesic_integrator/geodesic_integrator.hpp"  // GeodesicIntegrator
#include "../input_reader/input_reader.hpp"                // InputReader
#include "../utils/array.hpp"                              // Array
#include "../utils/exceptions.hpp"                         // BlacklightWarning

//--------------------------------------------------------------------------------------------------

// Radiation integrator constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
//   p_geodesic_integrator: pointer to object containing ray data
//   p_athena_reader_: pointer to object containing raw simulation data
RadiationIntegrator::RadiationIntegrator(const InputReader *p_input_reader,
    const GeodesicIntegrator *p_geodesic_integrator, const AthenaReader *p_athena_reader_)
  : p_athena_reader(p_athena_reader_)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();
  num_threads = p_input_reader->num_threads.value();

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

  // Copy checkpoint parameters
  if (model_type == ModelType::simulation)
  {
    checkpoint_sample_save = p_input_reader->checkpoint_sample_save.value();
    checkpoint_sample_load = p_input_reader->checkpoint_sample_load.value();
    if (checkpoint_sample_save and checkpoint_sample_load)
      throw BlacklightException("Cannot both save and load a sample checkpoint.");
    if (checkpoint_sample_save or checkpoint_sample_load)
      checkpoint_sample_file = p_input_reader->checkpoint_sample_file.value();
  }
  else
  {
    if (p_input_reader->checkpoint_sample_save.has_value()
        and p_input_reader->checkpoint_sample_save.value())
      BlacklightWarning("Ignoring checkpoint_sample_save selection.");
    if (p_input_reader->checkpoint_sample_load.has_value()
        and p_input_reader->checkpoint_sample_load.value())
      BlacklightWarning("Ignoring checkpoint_sample_load selection.");
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
  image_r = p_input_reader->image_r.value();
  image_resolution = p_input_reader->image_resolution.value();
  image_frequency = p_input_reader->image_frequency.value();
  if (model_type == ModelType::simulation)
    image_polarization = p_input_reader->image_polarization.value();
  else if (p_input_reader->image_polarization.has_value()
      and p_input_reader->image_polarization.value())
    BlacklightWarning("Ignoring image_polarization selection.");

  // Copy ray-tracing parameters
  ray_flat = p_input_reader->ray_flat.value();

  // Copy adaptive parameters
  adaptive_on = p_input_reader->adaptive_on.value();
  if (adaptive_on)
  {
    adaptive_block_size = p_input_reader->adaptive_block_size.value();
    if (adaptive_block_size <= 0)
      throw BlacklightException("Must have positive adaptive_block_size.");
    if (image_resolution % adaptive_block_size != 0)
      throw BlacklightException("Must have adaptive_block_size divide image_resolution.");
    adaptive_max_level = p_input_reader->adaptive_max_level.value();
    if (adaptive_max_level < 1)
      throw BlacklightException("Must have at least one allowed refinement level.");
    adaptive_val_frac = p_input_reader->adaptive_val_frac.value();
    if (adaptive_val_frac >= 0.0)
      adaptive_val_cut = p_input_reader->adaptive_val_cut.value();
    adaptive_abs_grad_frac = p_input_reader->adaptive_abs_grad_frac.value();
    if (adaptive_abs_grad_frac >= 0.0)
      adaptive_abs_grad_cut = p_input_reader->adaptive_abs_grad_cut.value();
    adaptive_rel_grad_frac = p_input_reader->adaptive_rel_grad_frac.value();
    if (adaptive_rel_grad_frac >= 0.0)
      adaptive_rel_grad_cut = p_input_reader->adaptive_rel_grad_cut.value();
    adaptive_abs_lapl_frac = p_input_reader->adaptive_abs_lapl_frac.value();
    if (adaptive_abs_lapl_frac >= 0.0)
      adaptive_abs_lapl_cut = p_input_reader->adaptive_abs_lapl_cut.value();
    adaptive_rel_lapl_frac = p_input_reader->adaptive_rel_lapl_frac.value();
    if (adaptive_rel_lapl_frac >= 0.0)
      adaptive_rel_lapl_cut = p_input_reader->adaptive_rel_lapl_cut.value();
  }

  // Copy camera data
  momentum_factor = p_geodesic_integrator->momentum_factor;
  for (int mu = 0; mu < 4; mu++)
  {
    camera_u_con[mu] = p_geodesic_integrator->u_con[mu];
    camera_u_cov[mu] = p_geodesic_integrator->u_cov[mu];
    camera_vert_con_c[mu] = p_geodesic_integrator->vert_con_c[mu];
  }
  camera_num_pix = p_geodesic_integrator->camera_num_pix;

  // Make shallow copies of camera arrays
  camera_pos = p_geodesic_integrator->camera_pos;
  camera_dir = p_geodesic_integrator->camera_dir;

  // Copy geodesic data
  geodesic_num_steps = p_geodesic_integrator->geodesic_num_steps;

  // Make shallow copies of geodesic arrays
  sample_flags = p_geodesic_integrator->sample_flags;
  sample_num = p_geodesic_integrator->sample_num;
  sample_pos = p_geodesic_integrator->sample_pos;
  sample_dir = p_geodesic_integrator->sample_dir;
  sample_len = p_geodesic_integrator->sample_len;

  // Copy adaptive geodesic data
  if (adaptive_on)
  {
    if (image_polarization)
    {
      camera_pos_adaptive = p_geodesic_integrator->camera_pos_adaptive;
      camera_dir_adaptive = p_geodesic_integrator->camera_dir_adaptive;
    }
    geodesic_num_steps_adaptive = p_geodesic_integrator->geodesic_num_steps_adaptive;
    sample_flags_adaptive = p_geodesic_integrator->sample_flags_adaptive;
    sample_num_adaptive = p_geodesic_integrator->sample_num_adaptive;
    sample_pos_adaptive = p_geodesic_integrator->sample_pos_adaptive;
    sample_dir_adaptive = p_geodesic_integrator->sample_dir_adaptive;
    sample_len_adaptive = p_geodesic_integrator->sample_len_adaptive;
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

  // Allocate space for calculating adaptive refinement
  if (adaptive_on)
  {
    linear_root_blocks = image_resolution / adaptive_block_size;
    block_num_pix = adaptive_block_size * adaptive_block_size;
    block_counts = new int[adaptive_max_level+1];
    block_counts[0] = linear_root_blocks * linear_root_blocks;
    refinement_flags = new Array<bool>[adaptive_max_level+1];
    refinement_flags[0].Allocate(block_counts[0]);
    image_adaptive = new Array<double>[adaptive_max_level+1];
    image_blocks = new Array<double>[num_threads];
    for (int thread = 0; thread < num_threads; thread++)
      if (model_type == ModelType::simulation and image_polarization)
        image_blocks[thread].Allocate(4, adaptive_block_size, adaptive_block_size);
      else
        image_blocks[thread].Allocate(adaptive_block_size, adaptive_block_size);
  }
}

//--------------------------------------------------------------------------------------------------

// Radiation integrator destructor
RadiationIntegrator::~RadiationIntegrator()
{
  if (adaptive_on)
  {
    for (int level = 0; level <= adaptive_max_level; level++)
    {
      refinement_flags[level].Deallocate();
      image_adaptive[level].Deallocate();
    }
    for (int thread = 0; thread < num_threads; thread++)
      image_blocks[thread].Deallocate();
    delete[] block_counts;
    delete[] refinement_flags;
    delete[] image_adaptive;
    delete[] image_blocks;
  }
}

//--------------------------------------------------------------------------------------------------

// Top-level function for processing raw data into image
// Inputs:
//   *p_time_sample: amount of time already taken for sampling
//   *p_time_integrate: amount of time already taken for integrating
// Outputs:
//   *p_time_sample: incremented by additional time taken for sampling
//   *p_time_integrate: incremented by additional time taken for integrating
//   returned value: flag indicating no additional geodesics need to be run for this snapshot
// Notes:
//   Assumes all data arrays have been set.
bool RadiationIntegrator::Integrate(double *p_time_sample, double *p_time_integrate)
{
  // Prepare timers
  double time_sample_start = 0.0;
  double time_sample_end = 0.0;
  double time_integrate_start = 0.0;
  double time_integrate_end = 0.0;

  // Sample simulation data
  if (model_type == ModelType::simulation)
  {
    time_sample_start = omp_get_wtime();
    if (first_time)
      ObtainGridData();
    if (adaptive_on and adaptive_current_level > 0)
      CalculateSimulationSampling();
    else if (first_time)
    {
      if (checkpoint_sample_load)
        LoadSampling();
      else
        CalculateSimulationSampling();
      if (checkpoint_sample_save)
        SaveSampling();
    }
    SampleSimulation();
    time_sample_end = omp_get_wtime();
  }

  // Integrate according to simulation data
  if (model_type == ModelType::simulation)
  {
    time_integrate_start = time_sample_end;
    CalculateSimulationCoefficients();
    if (image_polarization)
      IntegratePolarizedRadiation();
    else
      IntegrateUnpolarizedRadiation();
  }

  // Integrate according to formula
  if (model_type == ModelType::formula)
  {
    time_integrate_start = omp_get_wtime();
    CalculateFormulaCoefficients();
    IntegrateUnpolarizedRadiation();
  }

  // Check for adaptive refinement
  bool adaptive_complete = true;
  if (adaptive_on)
    adaptive_complete = CheckAdaptiveRefinement();
  if (adaptive_complete)
  {
    adaptive_num_levels = adaptive_current_level;
    adaptive_current_level = 0;
  }
  else
    adaptive_current_level++;
  time_integrate_end = omp_get_wtime();

  // Update first time flag
  first_time = false;

  // Calculate elapsed time
  *p_time_sample += time_sample_end - time_sample_start;
  *p_time_integrate += time_integrate_end - time_integrate_start;
  return adaptive_complete;
}
