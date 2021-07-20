// Blacklight radiation integrator

// C++ headers
#include <algorithm>  // max
#include <optional>   // optional

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
  image_resolution = p_input_reader->image_resolution.value();
  image_frequency = p_input_reader->image_frequency.value();
  if (model_type == ModelType::simulation)
    image_polarization = p_input_reader->image_polarization.value();
  else if (p_input_reader->image_polarization.has_value()
      and p_input_reader->image_polarization.value())
    BlacklightWarning("Ignoring image_polarization selection.");

  // Copy ray-tracing parameters
  ray_flat = p_input_reader->ray_flat.value();

  // Copy camera data
  momentum_factor = p_geodesic_integrator->momentum_factor;
  for (int mu = 0; mu < 4; mu++)
  {
    camera_ucon[mu] = p_geodesic_integrator->camera_ucon[mu];
    camera_ucov[mu] = p_geodesic_integrator->camera_ucov[mu];
    camera_up_con_c[mu] = p_geodesic_integrator->camera_up_con_c[mu];
  }

  // Make shallow copies of camera arrays
  image_position = p_geodesic_integrator->image_position;
  image_direction = p_geodesic_integrator->image_direction;

  // Copy geodesic data
  image_steps = p_geodesic_integrator->image_steps;

  // Make shallow copies of geodesic arrays
  sample_flags = p_geodesic_integrator->sample_flags;
  sample_num = p_geodesic_integrator->sample_num;
  sample_pos = p_geodesic_integrator->sample_pos;
  sample_dir = p_geodesic_integrator->sample_dir;
  sample_len = p_geodesic_integrator->sample_len;

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
}

//--------------------------------------------------------------------------------------------------

// Top-level function for processing raw data into image
// Inputs: (none)
// Output: (none)
// Notes:
//   Assumes all data arrays have been set.
void RadiationIntegrator::Integrate()
{
  // Calculate coefficients and integrate
  switch (model_type)
  {
    case ModelType::simulation:
    {
      if (first_time)
        CalculateSimulationSampling();
      SampleSimulation();
      CalculateSimulationCoefficients();
      if (image_polarization)
        IntegratePolarizedRadiation();
      else
        IntegrateUnpolarizedRadiation();
      break;
    }
    case ModelType::formula:
    {
      CalculateFormulaCoefficients();
      IntegrateUnpolarizedRadiation();
      break;
    }
  }

  // Update first time flag
  first_time = false;
  return;
}
