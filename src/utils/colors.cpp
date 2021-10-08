// Blacklight color utilities

// C++ headers
#include <cmath>  // pow

// Blacklight headers
#include "colors.hpp"

//--------------------------------------------------------------------------------------------------

// Function for converting from RGB to XYZ
// Inputs:
//   r, g, b: sRGB255 color coordinates (allowed to be fractional)
// Outputs:
//   *p_x, *p_y, *p_z: XYZ1 color coordinates under D65 illuminant
// Notes:
//   The conversion from sRGB1 to lRGB1 is the inverse of the (slightly non-monotonic) definition
//       sX = lX <= 0.0031308 ? 12.92 * lX : 1.055 * lX^(1/2.4) - 0.055, where the linear relation
//       is used to determine the transition in the other direction.
//   The matrix for converting lRGB1 to XYZ1 is the inverse of the definition
//       | 3.2406 -1.5372 -0.4986|
//       |-0.9689  1.8758  0.0415|
//       | 0.0557 -0.2040  1.0570|.
void RGBToXYZ(double r, double g, double b, double *p_x, double *p_y, double *p_z)
{
  double r1 = r / 255.0;
  double g1 = g / 255.0;
  double b1 = b / 255.0;
  double lr = r1 <= 0.040449936 ? r1 / 12.92 : std::pow((r1 + 0.055) / 1.055, 2.4);
  double lg = g1 <= 0.040449936 ? g1 / 12.92 : std::pow((g1 + 0.055) / 1.055, 2.4);
  double lb = b1 <= 0.040449936 ? b1 / 12.92 : std::pow((b1 + 0.055) / 1.055, 2.4);
  *p_x = 0.4123955889674142 * lr + 0.3575834307637148 * lg + 0.18049264738170154 * lb;
  *p_y = 0.21258623078559552 * lr + 0.715170303703411 * lg + 0.0722004986433362 * lb;
  *p_z = 0.019297215491746938 * lr + 0.11918386458084851 * lg + 0.9504971251315798 * lb;
  return;
}
