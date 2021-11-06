// Blacklight input reader - enum readers

// C++ headers
#include <string>  // string

// Blacklight headers
#include "input_reader.hpp"
#include "../blacklight.hpp"        // enums
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as ModelType enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid ModelType
// Notes:
//   Valid options:
//     "simulation": Athena++ output
//     "formula": parameterized formula from 2020 ApJ 897 148
ModelType InputReader::ReadModelType(const std::string &string)
{
  if (string == "simulation")
    return ModelType::simulation;
  else if (string == "formula")
    return ModelType::formula;
  else
    throw BlacklightException("Unknown string used for ModelType value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as OutputFormat enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid OutputFormat
// Notes:
//   Valid options:
//     "npz": NumPy .npz file with image and metadata arrays
//     "npy": NumPy .npy file with image array and minimal metadata
//     "raw": raw binary values of image array with no extra data
OutputFormat InputReader::ReadOutputFormat(const std::string &string)
{
  if (string == "npz")
    return OutputFormat::npz;
  else if (string == "npy")
    return OutputFormat::npy;
  else if (string == "raw")
    return OutputFormat::raw;
  else
    throw BlacklightException("Unknown string used for OutputFormat value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as Coordinates enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid Coordinates
// Notes:
//   Valid options:
//     "sph_ks": spherical Kerr-Schild
//     "cart_ks": Cartesian Kerr-Schild
Coordinates InputReader::ReadCoordinates(const std::string &string)
{
  if (string == "sph_ks")
    return Coordinates::sph_ks;
  else if (string == "cart_ks")
    return Coordinates::cart_ks;
  else
    throw BlacklightException("Unknown string used for Coordinates value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as Camera enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid Camera
// Notes:
//   Valid options:
//     "plane": rays originate from plane, parallel to central ray (which is perpendicular to
//         plane); appropriate for images as seen from infinity
//     "pinhole": rays originate from point, going in different directions; appropriate for camera
//         located near source
Camera InputReader::ReadCamera(const std::string &string)
{
  if (string == "plane")
    return Camera::plane;
  else if (string == "pinhole")
    return Camera::pinhole;
  else
    throw BlacklightException("Unknown string used for Camera value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as RayTerminate enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid RayTerminate
// Notes:
//   Valid options:
//     "photon": rays are terminated upon reaching the prograde equatorial photon orbit radius
//     "multiplicative": rays are terminated upon reaching the horizon radius times ray_factor
//         (dimensionless)
//     "additive": rays are terminated upon reaching the horizon radius plus ray_factor
//         (gravitational units)
RayTerminate InputReader::ReadRayTerminate(const std::string &string)
{
  if (string == "photon")
    return RayTerminate::photon;
  else if (string == "multiplicative")
    return RayTerminate::multiplicative;
  else if (string == "additive")
    return RayTerminate::additive;
  else
    throw BlacklightException("Unknown string used for RayTerminate value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as RayIntegrator enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid RayIntegrator
// Notes:
//   Valid options:
//     "dp": Dormand-Prince (5th-order adaptive Runge-Kutta)
//     "rk4": 4th-order Runge-Kutta
//     "rk2": 2nd-order Runge-Kutta
RayIntegrator InputReader::ReadRayIntegrator(const std::string &string)
{
  if (string == "dp")
    return RayIntegrator::dp;
  else if (string == "rk4")
    return RayIntegrator::rk4;
  else if (string == "rk2")
    return RayIntegrator::rk2;
  else
    throw BlacklightException("Unknown string used for RayIntegrator value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as FrequencySpacing enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid FrequencySpacing
// Notes:
//   Valid options:
//     "lin_freq": frequencies are linearly spaced from image_frequency_start to image_frequency_end
//     "lin_wave": wavelengths are linearly spaced from c / image_frequency_start to
//         c / image_frequency_end
//     "log": frequencies are logarithmically spaced from image_frequency_start to
//         image_frequency_end (and thus wavelengths are logarithmically spaced from
//         c / image_frequency_end to c / image_frequency_end)
FrequencySpacing InputReader::ReadFrequencySpacing(const std::string &string)
{
  if (string == "lin_freq")
    return FrequencySpacing::lin_freq;
  else if (string == "lin_wave")
    return FrequencySpacing::lin_wave;
  else if (string == "log")
    return FrequencySpacing::log;
  else
    throw BlacklightException("Unknown string used for FrequencySpacing value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as FrequencyNormalization enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid FrequencyNormalization
// Notes:
//   Valid options:
//     "camera": input image_frequency (or each frequency from image_frequency_min to
//         image_frequency_max) is taken to be the frequency as seen by the center of the camera,
//         accounting for its position and velocity; that is, the covariant time component of photon
//         momentum in the camera frame at the camera location is -image_frequency
//     "infinity": input image_frequency (or each frequency from image_frequency_min to
//         image_frequency_max) is taken to be the frequency of the light at the center of the
//         camera were it to be transported along geodesics to infinity and measured by an observer
//         at rest; that is, the covariant time component of photon momentum in the coordinate frame
//         is -image_frequency
FrequencyNormalization InputReader::ReadFrequencyNormalization(const std::string &string)
{
  if (string == "camera")
    return FrequencyNormalization::camera;
  else if (string == "infinity")
    return FrequencyNormalization::infinity;
  else
    throw BlacklightException("Unknown string used for FrequencyNormalization value.");
}

//--------------------------------------------------------------------------------------------------

// Function for interpreting strings as PlasmaModel enums
// Inputs:
//   string: string to be interpreted
// Outputs:
//   returned value: valid PlasmaModel
// Notes:
//   Valid options:
//     "ti_te_beta": spherical Kerr-Schild
//     "code_kappa": Cartesian Kerr-Schild
PlasmaModel InputReader::ReadPlasmaModel(const std::string &string)
{
  if (string == "ti_te_beta")
    return PlasmaModel::ti_te_beta;
  else if (string == "code_kappa")
    return PlasmaModel::code_kappa;
  else
    throw BlacklightException("Unknown string used for PlasmaModel value.");
}
