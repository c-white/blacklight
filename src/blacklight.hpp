// Blacklight main header

#ifndef BLACKLIGHT_H_
#define BLACKLIGHT_H_

// C++ headers
#include <complex>

// Mathematical constants
namespace Math
{
  constexpr double pi = 3.141592653589793;
  constexpr double sqrt2 = 1.4142135623730951;
  constexpr std::complex<double> i(0.0, 1.0);
}

// Physical constants
namespace Physics
{
  constexpr double c = 2.99792458e10;
  constexpr double h = 6.62607015e-27;
  constexpr double k_b = 1.380649e-16;
  constexpr double m_p = 1.67262192369e-24;
  constexpr double m_e = 9.1093837015e-28;
  constexpr double e = 4.80320425e-10;
  constexpr double gg_msun = 1.32712440018e26;
}

// Unscoped enumerations
namespace CellValues
{
  enum : int {rho, n_e, p_gas, theta_e, bb, sigma, beta_inv, num_cell_values};
}

// Scoped enumerations
enum struct ModelType {simulation, formula};
enum struct OutputFormat {raw, npy, npz};
enum struct Coordinates {sph_ks, cart_ks};
enum struct Camera {plane, pinhole};
enum struct RayTerminate {photon, multiplicative, additive};
enum struct FrequencySpacing {lin_freq, lin_wave, log};
enum struct FrequencyNormalization {camera, infinity};
enum struct RenderType {fill, thresh, rise, fall};
enum struct PlasmaModel {ti_te_beta, code_kappa};

#endif
