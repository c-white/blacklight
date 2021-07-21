// Blacklight geodesic integrator

// C++ headers
#include <cmath>     // abs, acos, cos, sqrt
#include <optional>  // optional

// Library headers
#include <omp.h>  // omp_get_wtime

// Blacklight headers
#include "geodesic_integrator.hpp"
#include "../input_reader/input_reader.hpp"  // InputReader

//--------------------------------------------------------------------------------------------------

// Geodesic integrator constructor
// Inputs:
//   p_input_reader: pointer to object containing input parameters
GeodesicIntegrator::GeodesicIntegrator(const InputReader *p_input_reader)
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

  // Calculate number of pixels
  camera_num_pix = image_resolution * image_resolution;
}

//--------------------------------------------------------------------------------------------------

// Top-level function for integrating geodesics
// Inputs: (none)
// Output:
//   returned value: execution time in seconds
// Notes:
//   Assumes all data arrays have been set.
double GeodesicIntegrator::Integrate()
{
  double time_start = omp_get_wtime();
  InitializeCamera();
  InitializeGeodesics();
  IntegrateGeodesics();
  ReverseGeodesics();
  return omp_get_wtime() - time_start;
}
