// Blacklight main header

#ifndef BLACKLIGHT_H_
#define BLACKLIGHT_H_

// Mathematical constants
namespace math
{
  constexpr double pi = 3.141592653589793;
}

// Physical constants
namespace physics
{
  constexpr double c = 2.99792458e10;
  constexpr double gg_msun = 1.32712440018e26;
  constexpr double h = 6.62607015e-27;
  constexpr double k_b = 1.380649e-16;
  constexpr double e = 4.80320425e-10;
  constexpr double m_p = 1.67262192369e-24;
  constexpr double m_e = 9.1093837015e-28;
}

// Enumerations
enum Coordinates {sph_ks, cart_ks};

#endif
