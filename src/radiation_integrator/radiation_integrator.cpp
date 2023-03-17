// Blacklight radiation integrator

// C++ headers
#include <algorithm>  // max
#include <optional>   // optional

// Library headers
#include <omp.h>  // omp_get_wtime

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"                               // Physics, enums
#include "../geodesic_integrator/geodesic_integrator.hpp"  // GeodesicIntegrator
#include "../input_reader/input_reader.hpp"                // InputReader
#include "../simulation_reader/simulation_reader.hpp"      // SimulationReader
#include "../utils/array.hpp"                              // Array
#include "../utils/exceptions.hpp"                         // BlacklightException, BlacklightWarning

//--------------------------------------------------------------------------------------------------

// Radiation integrator constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
//   p_geodesic_integrator: pointer to object containing ray data
//   p_simulation_reader_: pointer to object containing raw simulation data
RadiationIntegrator::RadiationIntegrator(const InputReader *p_input_reader,
    const GeodesicIntegrator *p_geodesic_integrator, const SimulationReader *p_simulation_reader_)
  : p_simulation_reader(p_simulation_reader_)
{
  // Copy general input data
  model_type = p_input_reader->model_type.value();
  num_threads = p_input_reader->num_threads.value();

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

  // Copy simulation parameters
  if (model_type == ModelType::simulation)
  {
    simulation_format = p_input_reader->simulation_format.value();
    simulation_coord = p_input_reader->simulation_coord.value();
    simulation_m_msun = p_input_reader->simulation_m_msun.value();
    simulation_rho_cgs = p_input_reader->simulation_rho_cgs.value();
    simulation_interp = p_input_reader->simulation_interp.value();
    if ((simulation_format == SimulationFormat::athena
        or simulation_format == SimulationFormat::athenak) and simulation_interp)
      simulation_block_interp = p_input_reader->simulation_block_interp.value();
    else if (p_input_reader->simulation_block_interp.has_value())
      BlacklightWarning("Ignoring simulation_block_interp selection.");
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

  // Copy camera parameters
  camera_r = p_input_reader->camera_r.value();
  camera_width = p_input_reader->camera_width.value();
  camera_resolution = p_input_reader->camera_resolution.value();
  if (camera_resolution <= 0)
    throw BlacklightException("Must have positive camera_resolution.");

  // Copy ray-tracing parameters
  ray_flat = p_input_reader->ray_flat.value();

  // Copy image parameters
  image_light = p_input_reader->image_light.value();
  image_num_frequencies = p_input_reader->image_num_frequencies.value();
  if (image_light)
  {
    if (model_type == ModelType::simulation)
      image_polarization = p_input_reader->image_polarization.value();
    else if (p_input_reader->image_polarization.has_value()
        and p_input_reader->image_polarization.value())
      BlacklightWarning("Ignoring image_polarization selection.");
    if (image_polarization)
      image_rotation_split = p_input_reader->image_rotation_split.value();
  }
  else if (p_input_reader->image_polarization.has_value()
      and p_input_reader->image_polarization.value())
    BlacklightWarning("Ignoring image_polarization selection.");
  image_time = p_input_reader->image_time.value();
  image_length = p_input_reader->image_length.value();
  image_lambda = p_input_reader->image_lambda.value();
  image_emission = p_input_reader->image_emission.value();
  image_tau = p_input_reader->image_tau.value();
  if (model_type == ModelType::simulation)
  {
    image_lambda_ave = p_input_reader->image_lambda_ave.value();
    image_emission_ave = p_input_reader->image_emission_ave.value();
    image_tau_int = p_input_reader->image_tau_int.value();
  }
  else
  {
    if (p_input_reader->image_lambda_ave.has_value() and p_input_reader->image_lambda_ave.value())
      BlacklightWarning("Ignoring image_lambda_ave selection.");
    image_lambda_ave = false;
    if (p_input_reader->image_emission_ave.has_value()
        and p_input_reader->image_emission_ave.value())
      BlacklightWarning("Ignoring image_emission_ave selection.");
    image_emission_ave = false;
    if (p_input_reader->image_tau_int.has_value() and p_input_reader->image_tau_int.value())
      BlacklightWarning("Ignoring image_tau_int selection.");
    image_tau_int = false;
  }
  image_crossings = p_input_reader->image_crossings.value();

  // Copy rendering parameters
  if (model_type == ModelType::simulation)
    render_num_images = p_input_reader->render_num_images.value();
  else
  {
    if (p_input_reader->render_num_images.has_value()
        and p_input_reader->render_num_images.value() > 0)
      BlacklightWarning("Ignoring request for rendering.");
    render_num_images = 0;
  }
  if (render_num_images > 0)
  {
    render_num_features = new int[render_num_images];
    render_quantities = new int *[render_num_images]();
    render_types = new RenderType *[render_num_images]();
    render_min_vals = new double *[render_num_images]();
    render_max_vals = new double *[render_num_images]();
    render_thresh_vals = new double *[render_num_images]();
    render_tau_scales = new double *[render_num_images]();
    render_opacities = new double *[render_num_images]();
    render_x_vals = new double *[render_num_images]();
    render_y_vals = new double *[render_num_images]();
    render_z_vals = new double *[render_num_images]();
    for (int n_i = 0; n_i < render_num_images; n_i++)
    {
      int num_features = p_input_reader->render_num_features[n_i].value();
      if (num_features <= 0)
        throw BlacklightException("Must have positive number of features for each rendered image.");
      render_num_features[n_i] = num_features;
      render_quantities[n_i] = new int[num_features];
      render_types[n_i] = new RenderType[num_features];
      render_min_vals[n_i] = new double[num_features];
      render_max_vals[n_i] = new double[num_features];
      render_thresh_vals[n_i] = new double[num_features];
      render_tau_scales[n_i] = new double[num_features];
      render_opacities[n_i] = new double[num_features];
      render_x_vals[n_i] = new double[num_features];
      render_y_vals[n_i] = new double[num_features];
      render_z_vals[n_i] = new double[num_features];
      for (int n_f = 0; n_f < num_features; n_f++)
      {
        render_quantities[n_i][n_f] = p_input_reader->render_quantities[n_i][n_f].value();
        render_types[n_i][n_f] = p_input_reader->render_types[n_i][n_f].value();
        if (render_types[n_i][n_f] == RenderType::fill)
        {
          render_min_vals[n_i][n_f] = p_input_reader->render_min_vals[n_i][n_f].value();
          render_max_vals[n_i][n_f] = p_input_reader->render_max_vals[n_i][n_f].value();
          render_tau_scales[n_i][n_f] = p_input_reader->render_tau_scales[n_i][n_f].value();
        }
        if (render_types[n_i][n_f] == RenderType::thresh
            or render_types[n_i][n_f] == RenderType::rise
            or render_types[n_i][n_f] == RenderType::fall)
        {
          render_thresh_vals[n_i][n_f] = p_input_reader->render_thresh_vals[n_i][n_f].value();
          render_opacities[n_i][n_f] = p_input_reader->render_opacities[n_i][n_f].value();
        }
        render_x_vals[n_i][n_f] = p_input_reader->render_x_vals[n_i][n_f].value();
        render_y_vals[n_i][n_f] = p_input_reader->render_y_vals[n_i][n_f].value();
        render_z_vals[n_i][n_f] = p_input_reader->render_z_vals[n_i][n_f].value();
      }
    }
  }
  if (not (image_light or image_time or image_length or image_lambda or image_emission or image_tau
      or image_lambda_ave or image_emission_ave or image_tau_int or image_crossings
      or render_num_images > 0))
    throw BlacklightException("No image or rendering selected.");

  // Copy slow-light parameters
  if (model_type == ModelType::simulation)
  {
    slow_light_on = p_input_reader->slow_light_on.value();
    if (slow_light_on)
    {
      if (checkpoint_sample_save or checkpoint_sample_load)
        throw BlacklightException("Cannot use sample checkpoints with slow light.");
      slow_interp = p_input_reader->slow_interp.value();
      slow_chunk_size = p_input_reader->slow_chunk_size.value();
      slow_t_start = p_input_reader->slow_t_start.value();
      slow_dt = p_input_reader->slow_dt.value();
    }
  }

  // Copy adaptive parameters
  adaptive_max_level = p_input_reader->adaptive_max_level.value();
  if (adaptive_max_level > 0)
  {
    if (not image_light)
      throw BlacklightException("Adaptive ray tracing requires image_light.");
    adaptive_block_size = p_geodesic_integrator->adaptive_block_size;
    if (adaptive_block_size <= 0)
      throw BlacklightException("Must have positive adaptive_block_size.");
    if (camera_resolution % adaptive_block_size != 0)
      throw BlacklightException("Must have adaptive_block_size divide camera_resolution.");
    if (image_num_frequencies > 1)
    {
      adaptive_frequency_num = p_input_reader->adaptive_frequency_num.value() - 1;
      if (adaptive_frequency_num < 0 or adaptive_frequency_num >= image_num_frequencies)
        throw BlacklightException(
            "Must choose adaptive_frequency_num from 1 to image_num_frequencies.");
    }
    else
      adaptive_frequency_num = 0;
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
    adaptive_num_regions = p_input_reader->adaptive_num_regions.value();
    if (adaptive_num_regions > 0)
    {
      adaptive_region_levels = new int[adaptive_num_regions];
      adaptive_region_x_min_vals = new double[adaptive_num_regions];
      adaptive_region_x_max_vals = new double[adaptive_num_regions];
      adaptive_region_y_min_vals = new double[adaptive_num_regions];
      adaptive_region_y_max_vals = new double[adaptive_num_regions];
      for (int n_r = 0; n_r < adaptive_num_regions; n_r++)
      {
        adaptive_region_levels[n_r] = p_input_reader->adaptive_region_levels[n_r].value();
        adaptive_region_x_min_vals[n_r] = p_input_reader->adaptive_region_x_min_vals[n_r].value();
        adaptive_region_x_max_vals[n_r] = p_input_reader->adaptive_region_x_max_vals[n_r].value();
        adaptive_region_y_min_vals[n_r] = p_input_reader->adaptive_region_y_min_vals[n_r].value();
        adaptive_region_y_max_vals[n_r] = p_input_reader->adaptive_region_y_max_vals[n_r].value();
      }
    }
  }

  // Copy plasma parameters
  if (model_type == ModelType::simulation)
  {
    plasma_mu = p_input_reader->plasma_mu.value();
    plasma_ne_ni = p_input_reader->plasma_ne_ni.value();
    plasma_model = p_input_reader->plasma_model.value();
    if (plasma_model == PlasmaModel::ti_te_beta)
    {
      plasma_rat_low = p_input_reader->plasma_rat_low.value();
      plasma_rat_high = p_input_reader->plasma_rat_high.value();
    }
    plasma_power_frac = p_input_reader->plasma_power_frac.value();
    if (plasma_power_frac < 0.0 or plasma_power_frac > 1.0)
      BlacklightWarning("Fraction of power-law electrons outside [0, 1].");
    if (plasma_power_frac != 0.0)
    {
      plasma_p = p_input_reader->plasma_p.value();
      plasma_gamma_min = p_input_reader->plasma_gamma_min.value();
      plasma_gamma_max = p_input_reader->plasma_gamma_max.value();
    }
    plasma_kappa_frac = p_input_reader->plasma_kappa_frac.value();
    if (plasma_kappa_frac < 0.0 or plasma_kappa_frac > 1.0)
      BlacklightWarning("Fraction of kappa-distribution electrons outside [0, 1].");
    if (plasma_kappa_frac != 0.0)
    {
      plasma_kappa = p_input_reader->plasma_kappa.value();
      if (image_light and image_polarization)
      {
        if (plasma_kappa < 3.5 or plasma_kappa > 5.0)
          throw BlacklightException("Polarized transport only supports kappa in [3.5, 5].");
        else if (plasma_kappa != 3.5 and plasma_kappa != 4.0 and plasma_kappa != 4.5
            and plasma_kappa != 5.0)
          BlacklightWarning("Polarized transport will interpolate formulas based on kappa.");
      }
      plasma_w = p_input_reader->plasma_w.value();
    }
    plasma_thermal_frac = 1.0 - (plasma_power_frac + plasma_kappa_frac);
    if (plasma_thermal_frac < 0.0 or plasma_thermal_frac > 1.0)
      BlacklightWarning("Fraction of thermal electrons outside [0, 1].");
  }

  // Copy cut parameters
  if (model_type == ModelType::simulation)
  {
    cut_rho_min = p_input_reader->cut_rho_min.value();
    cut_rho_max = p_input_reader->cut_rho_max.value();
    cut_n_e_min = p_input_reader->cut_n_e_min.value();
    cut_n_e_max = p_input_reader->cut_n_e_max.value();
    cut_p_gas_min = p_input_reader->cut_p_gas_min.value();
    cut_p_gas_max = p_input_reader->cut_p_gas_max.value();
    cut_theta_e_min = p_input_reader->cut_theta_e_min.value();
    cut_theta_e_max = p_input_reader->cut_theta_e_max.value();
    cut_b_min = p_input_reader->cut_b_min.value();
    cut_b_max = p_input_reader->cut_b_max.value();
    cut_sigma_min = p_input_reader->cut_sigma_min.value();
    cut_sigma_max = p_input_reader->cut_sigma_max.value();
    cut_beta_inverse_min = p_input_reader->cut_beta_inverse_min.value();
    cut_beta_inverse_max = p_input_reader->cut_beta_inverse_max.value();
  }
  cut_omit_near = p_input_reader->cut_omit_near.value();
  cut_omit_far = p_input_reader->cut_omit_far.value();
  cut_omit_in = p_input_reader->cut_omit_in.value();
  cut_omit_out = p_input_reader->cut_omit_out.value();
  cut_midplane_theta = p_input_reader->cut_midplane_theta.value();
  cut_midplane_z = p_input_reader->cut_midplane_z.value();
  cut_plane = p_input_reader->cut_plane.value();
  if (cut_plane)
  {
    cut_plane_origin_x = p_input_reader->cut_plane_origin_x.value();
    cut_plane_origin_y = p_input_reader->cut_plane_origin_y.value();
    cut_plane_origin_z = p_input_reader->cut_plane_origin_z.value();
    cut_plane_normal_x = p_input_reader->cut_plane_normal_x.value();
    cut_plane_normal_y = p_input_reader->cut_plane_normal_y.value();
    cut_plane_normal_z = p_input_reader->cut_plane_normal_z.value();
  }

  // Copy fallback parameters
  fallback_nan = p_input_reader->fallback_nan.value();
  if (model_type == ModelType::simulation and not fallback_nan)
  {
    fallback_rho = p_input_reader->fallback_rho.value();
    fallback_pgas = p_input_reader->fallback_pgas.value();
    if (plasma_model == PlasmaModel::code_kappa)
      fallback_kappa = p_input_reader->fallback_kappa.value();
  }

  // Copy camera data
  for (int mu = 0; mu < 4; mu++)
  {
    camera_x[mu] = p_geodesic_integrator->cam_x[mu];
    camera_u_con[mu] = p_geodesic_integrator->u_con[mu];
    camera_u_cov[mu] = p_geodesic_integrator->u_cov[mu];
    camera_vert_con_c[mu] = p_geodesic_integrator->vert_con_c[mu];
  }
  camera_num_pix = p_geodesic_integrator->camera_num_pix;

  // Make shallow copies of camera arrays
  camera_loc = p_geodesic_integrator->camera_loc;
  camera_pos = p_geodesic_integrator->camera_pos;
  camera_dir = p_geodesic_integrator->camera_dir;

  // Make shallow copies of image frequency arrays
  image_frequencies = p_geodesic_integrator->image_frequencies;
  momentum_factors = p_geodesic_integrator->momentum_factors;

  // Make shallow copy of geodesic data
  geodesic_num_steps = p_geodesic_integrator->geodesic_num_steps;

  // Make shallow copies of geodesic arrays
  sample_flags = p_geodesic_integrator->sample_flags;
  sample_num = p_geodesic_integrator->sample_num;
  sample_pos = p_geodesic_integrator->sample_pos;
  sample_dir = p_geodesic_integrator->sample_dir;
  sample_len = p_geodesic_integrator->sample_len;

  // Allocate space for sample data
  sample_inds = new Array<int>[adaptive_max_level+1];
  sample_fracs = new Array<double>[adaptive_max_level+1];
  sample_nan = new Array<bool>[adaptive_max_level+1];
  sample_cut = new Array<bool>[adaptive_max_level+1];
  sample_fallback = new Array<bool>[adaptive_max_level+1];
  sample_rho = new Array<float>[adaptive_max_level+1];
  sample_pgas = new Array<float>[adaptive_max_level+1];
  sample_kappa = new Array<float>[adaptive_max_level+1];
  sample_uu1 = new Array<float>[adaptive_max_level+1];
  sample_uu2 = new Array<float>[adaptive_max_level+1];
  sample_uu3 = new Array<float>[adaptive_max_level+1];
  sample_bb1 = new Array<float>[adaptive_max_level+1];
  sample_bb2 = new Array<float>[adaptive_max_level+1];
  sample_bb3 = new Array<float>[adaptive_max_level+1];

  // Allocate space for coefficient data
  j_i = new Array<double>[adaptive_max_level+1];
  j_q = new Array<double>[adaptive_max_level+1];
  j_v = new Array<double>[adaptive_max_level+1];
  alpha_i = new Array<double>[adaptive_max_level+1];
  alpha_q = new Array<double>[adaptive_max_level+1];
  alpha_v = new Array<double>[adaptive_max_level+1];
  rho_q = new Array<double>[adaptive_max_level+1];
  rho_v = new Array<double>[adaptive_max_level+1];
  cell_values = new Array<double>[adaptive_max_level+1];

  // Copy slow-light extrapolation tolerance
  if (slow_light_on)
    extrapolation_tolerance = p_simulation_reader->extrapolation_tolerance;

  // Set and calculate geometry data
  if (model_type == ModelType::simulation)
  {
    bh_m = 1.0;
    bh_a = p_input_reader->simulation_a.value();
    mass_msun = simulation_m_msun;
  }
  else if (model_type == ModelType::formula)
  {
    bh_m = 1.0;
    bh_a = p_input_reader->formula_spin.value();
    mass_msun = formula_mass * Physics::c * Physics::c / Physics::gg_msun;
  }

  // Allocate space for image data
  image = new Array<double>[adaptive_max_level+1];

  // Calculate number of simultaneous images and their offsets
  if (image_light)
  {
    image_num_quantities += image_num_frequencies
        * (model_type == ModelType::simulation and image_polarization ? 4 : 1);
    image_offset_time = image_num_quantities;
    image_offset_length = image_num_quantities;
    image_offset_lambda = image_num_quantities;
    image_offset_emission = image_num_quantities;
    image_offset_tau = image_num_quantities;
    image_offset_lambda_ave = image_num_quantities;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_time)
  {
    image_num_quantities++;
    image_offset_length = image_num_quantities;
    image_offset_lambda = image_num_quantities;
    image_offset_emission = image_num_quantities;
    image_offset_tau = image_num_quantities;
    image_offset_lambda_ave = image_num_quantities;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_length)
  {
    image_num_quantities++;
    image_offset_lambda = image_num_quantities;
    image_offset_emission = image_num_quantities;
    image_offset_tau = image_num_quantities;
    image_offset_lambda_ave = image_num_quantities;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_lambda)
  {
    image_num_quantities += image_num_frequencies;
    image_offset_emission = image_num_quantities;
    image_offset_tau = image_num_quantities;
    image_offset_lambda_ave = image_num_quantities;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_emission)
  {
    image_num_quantities += image_num_frequencies;
    image_offset_tau = image_num_quantities;
    image_offset_lambda_ave = image_num_quantities;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_tau)
  {
    image_num_quantities += image_num_frequencies;
    image_offset_lambda_ave = image_num_quantities;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_lambda_ave)
  {
    image_num_quantities += image_num_frequencies * CellValues::num_cell_values;
    image_offset_emission_ave = image_num_quantities;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_emission_ave)
  {
    image_num_quantities += image_num_frequencies * CellValues::num_cell_values;
    image_offset_tau_int = image_num_quantities;
    image_offset_crossings = image_num_quantities;
  }
  if (image_tau_int)
  {
    image_num_quantities += image_num_frequencies * CellValues::num_cell_values;
    image_offset_crossings = image_num_quantities;
  }
  if (image_crossings)
    image_num_quantities++;

  // Allocate space for rendering data
  render = new Array<double>[adaptive_max_level+1];

  // Allocate space for calculating adaptive refinement
  if (adaptive_max_level > 0)
  {
    linear_root_blocks = camera_resolution / adaptive_block_size;
    block_num_pix = adaptive_block_size * adaptive_block_size;
    block_counts = new int[adaptive_max_level+1];
    block_counts[0] = linear_root_blocks * linear_root_blocks;
    refinement_flags = new Array<bool>[adaptive_max_level+1];
    refinement_flags[0].Allocate(block_counts[0]);
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
  // Free memory - rendering input
  for (int n_i = 0; n_i < render_num_images; n_i++)
  {
    delete[] render_quantities[n_i];
    delete[] render_types[n_i];
    delete[] render_min_vals[n_i];
    delete[] render_max_vals[n_i];
    delete[] render_thresh_vals[n_i];
    delete[] render_tau_scales[n_i];
    delete[] render_opacities[n_i];
    delete[] render_x_vals[n_i];
    delete[] render_y_vals[n_i];
    delete[] render_z_vals[n_i];
  }
  delete[] render_num_features;
  delete[] render_quantities;
  delete[] render_types;
  delete[] render_min_vals;
  delete[] render_max_vals;
  delete[] render_thresh_vals;
  delete[] render_tau_scales;
  delete[] render_opacities;
  delete[] render_x_vals;
  delete[] render_y_vals;
  delete[] render_z_vals;

  // Free memory - adaptive input
  delete[] adaptive_region_levels;
  delete[] adaptive_region_x_min_vals;
  delete[] adaptive_region_x_max_vals;
  delete[] adaptive_region_y_min_vals;
  delete[] adaptive_region_y_max_vals;

  // Free memory - sample data
  for (int level = 0; level <= adaptive_max_level; level++)
  {
    sample_inds[level].Deallocate();
    sample_fracs[level].Deallocate();
    sample_nan[level].Deallocate();
    sample_cut[level].Deallocate();
    sample_fallback[level].Deallocate();
    sample_rho[level].Deallocate();
    sample_pgas[level].Deallocate();
    sample_kappa[level].Deallocate();
    sample_uu1[level].Deallocate();
    sample_uu2[level].Deallocate();
    sample_uu3[level].Deallocate();
    sample_bb1[level].Deallocate();
    sample_bb2[level].Deallocate();
    sample_bb3[level].Deallocate();
  }
  delete[] sample_inds;
  delete[] sample_fracs;
  delete[] sample_nan;
  delete[] sample_cut;
  delete[] sample_fallback;
  delete[] sample_rho;
  delete[] sample_pgas;
  delete[] sample_kappa;
  delete[] sample_uu1;
  delete[] sample_uu2;
  delete[] sample_uu3;
  delete[] sample_bb1;
  delete[] sample_bb2;
  delete[] sample_bb3;

  // Free memory - coefficient data
  for (int level = 0; level <= adaptive_max_level; level++)
  {
    j_i[level].Deallocate();
    j_q[level].Deallocate();
    j_v[level].Deallocate();
    alpha_i[level].Deallocate();
    alpha_q[level].Deallocate();
    alpha_v[level].Deallocate();
    rho_q[level].Deallocate();
    rho_v[level].Deallocate();
    cell_values[level].Deallocate();
  }
  delete[] j_i;
  delete[] j_q;
  delete[] j_v;
  delete[] alpha_i;
  delete[] alpha_q;
  delete[] alpha_v;
  delete[] rho_q;
  delete[] rho_v;
  delete[] cell_values;

  // Free memory - image data
  for (int level = 0; level <= adaptive_max_level; level++)
    image[level].Deallocate();
  delete[] image;

  // Free memory - rendering data
  for (int level = 0; level <= adaptive_max_level; level++)
    render[level].Deallocate();
  delete[] render;

  // Free memory - adaptive data
  if (adaptive_max_level > 0)
  {
    for (int level = 0; level <= adaptive_max_level; level++)
      refinement_flags[level].Deallocate();
    for (int thread = 0; thread < num_threads; thread++)
      image_blocks[thread].Deallocate();
    delete[] block_counts;
    delete[] refinement_flags;
    delete[] image_blocks;
  }
}

//--------------------------------------------------------------------------------------------------

// Top-level function for processing raw data into image and/or rendering
// Inputs:
//   snapshot: index (starting at 0) of which snapshot is about to be processed
//   *p_time_sample: amount of time already taken for sampling
//   *p_time_image: amount of time already taken for integrating image and evaluating refinement
//   *p_time_render: amount of time already taken for rendering
// Outputs:
//   *p_time_sample: incremented by additional time taken for sampling
//   *p_time_image: incremented by additional time taken for integrating image and evaluating
//       refinement
//   *p_time_render: incremented by additional time taken for rendering
//   returned value: flag indicating no additional geodesics need to be run for this snapshot
// Notes:
//   Assumes all data arrays have been set.
bool RadiationIntegrator::Integrate(int snapshot, double *p_time_sample, double *p_time_image,
    double *p_time_render)
{
  // Prepare timers
  double time_sample_start = 0.0;
  double time_sample_end = 0.0;
  double time_image_start = 0.0;
  double time_image_end = 0.0;
  double time_render_start = 0.0;
  double time_render_end = 0.0;
  double time_refine_start = 0.0;
  double time_refine_end = 0.0;

  // Sample simulation data
  if (model_type == ModelType::simulation)
  {
    time_sample_start = omp_get_wtime();
    if (first_time)
      ObtainGridData();
    if (adaptive_level > 0)
      CalculateSimulationSampling(snapshot);
    else if (first_time)
    {
      if (checkpoint_sample_load)
        LoadSampling();
      else
        CalculateSimulationSampling(snapshot);
      if (checkpoint_sample_save)
        SaveSampling();
    }
    else if (slow_light_on)
      CalculateSimulationSampling(snapshot);
    SampleSimulation();
    time_sample_end = omp_get_wtime();
  }

  // Integrate according to simulation data
  if (model_type == ModelType::simulation)
  {
    time_image_start = time_sample_end;
    CalculateSimulationCoefficients();
    if (image_light and image_polarization)
      IntegratePolarizedRadiation();
    else if (image_light or image_time or image_length or image_lambda or image_emission
        or image_tau or image_lambda_ave or image_emission_ave or image_tau_int or image_crossings)
      IntegrateUnpolarizedRadiation();
    time_image_end = omp_get_wtime();
    if (render_num_images > 0)
    {
      time_render_start = time_image_end;
      Render();
      time_render_end = omp_get_wtime();
    }
  }

  // Integrate according to formula
  if (model_type == ModelType::formula)
  {
    time_image_start = omp_get_wtime();
    CalculateFormulaCoefficients();
    IntegrateUnpolarizedRadiation();
    time_image_end = omp_get_wtime();
  }

  // Check for adaptive refinement
  time_refine_start = omp_get_wtime();
  bool adaptive_complete = true;
  if (adaptive_max_level > 0)
    adaptive_complete = CheckAdaptiveRefinement();
  if (adaptive_complete)
  {
    adaptive_num_levels = adaptive_level;
    adaptive_level = 0;
  }
  else
    adaptive_level++;
  time_refine_end = omp_get_wtime();

  // Update first time flag
  first_time = false;

  // Calculate elapsed time
  *p_time_sample += time_sample_end - time_sample_start;
  *p_time_image += time_image_end - time_image_start + time_refine_end - time_refine_start;
  *p_time_render += time_render_end - time_render_start;
  return adaptive_complete;
}
