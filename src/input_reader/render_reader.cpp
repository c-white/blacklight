// Blacklight input reader - render reader

// C++ headers
#include <cstddef>   // size_t
#include <optional>  // optional
#include <sstream>   // ostringstream
#include <string>    // stod, stoi, string

// Blacklight headers
#include "input_reader.hpp"
#include "../blacklight.hpp"        // enums
#include "../utils/colors.hpp"      // RGBToXYZ
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Function for parsing rendering options
// Inputs:
//   key: input key as a string without leading "render_"
//   val: input value as a string
// Outputs: (none)
// Notes:
//   Values are only recorded if space is allocated for them, so the keys should occur in a sensible
//       order in the input file.
//   Variables indexed beyond what is allocated (e.g. render_2_num_features, when
//       render_num_features = 1) are silently ignored.
void InputReader::ReadRender(const std::string &key, const std::string &val)
{
  // Read total number of render images
  if (key == "num_images")
  {
    render_num_images = std::stoi(val);
    if (render_num_images.value() > 0)
    {
      render_num_features = new std::optional<int>[render_num_images.value()];
      render_types = new std::optional<RenderType> *[render_num_images.value()]();
      render_quantities = new std::optional<int> *[render_num_images.value()]();
      render_min_vals = new std::optional<double> *[render_num_images.value()]();
      render_max_vals = new std::optional<double> *[render_num_images.value()]();
      render_thresh_vals = new std::optional<double> *[render_num_images.value()]();
      render_r_vals = new std::optional<double> *[render_num_images.value()]();
      render_tau_scales = new std::optional<double> *[render_num_images.value()]();
      render_opacities = new std::optional<double> *[render_num_images.value()]();
      render_x_vals = new std::optional<double> *[render_num_images.value()]();
      render_y_vals = new std::optional<double> *[render_num_images.value()]();
      render_z_vals = new std::optional<double> *[render_num_images.value()]();
      render_lx_vals = new std::optional<double> *[render_num_images.value()]();
      render_ly_vals = new std::optional<double> *[render_num_images.value()]();
      render_lz_vals = new std::optional<double> *[render_num_images.value()]();
      render_stream_files = new std::optional<std::string> *[render_num_images.value()]();
    }
  }

  // Read number of features in a particular image
  else if (key.size() >= 14 and key.compare(key.size() - 13, key.npos, "_num_features") == 0)
  {
    int image_num = std::stoi(key.substr(0, key.size() - 13)) - 1;
    if (image_num >= render_num_images.value())
      return;
    int num_features = std::stoi(val);
    render_num_features[image_num] = num_features;
    render_types[image_num] = new std::optional<RenderType>[num_features];
    render_quantities[image_num] = new std::optional<int>[num_features];
    render_min_vals[image_num] = new std::optional<double>[num_features];
    render_max_vals[image_num] = new std::optional<double>[num_features];
    render_thresh_vals[image_num] = new std::optional<double>[num_features];
    render_r_vals[image_num] = new std::optional<double>[num_features];
    render_tau_scales[image_num] = new std::optional<double>[num_features];
    render_opacities[image_num] = new std::optional<double>[num_features];
    render_x_vals[image_num] = new std::optional<double>[num_features];
    render_y_vals[image_num] = new std::optional<double>[num_features];
    render_z_vals[image_num] = new std::optional<double>[num_features];
    render_lx_vals[image_num] = new std::optional<double>[num_features];
    render_ly_vals[image_num] = new std::optional<double>[num_features];
    render_lz_vals[image_num] = new std::optional<double>[num_features];
    render_stream_files[image_num] = new std::optional<std::string>[num_features];
  }

  // Read type
  else if (key.size() >= 8 and key.compare(key.size() - 5, key.npos, "_type") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    if (val == "fill")
      render_types[image_num][feature_num] = RenderType::fill;
    else if (val == "thresh")
      render_types[image_num][feature_num] = RenderType::thresh;
    else if (val == "rise")
      render_types[image_num][feature_num] = RenderType::rise;
    else if (val == "fall")
      render_types[image_num][feature_num] = RenderType::fall;
    else if (val == "line")
      render_types[image_num][feature_num] = RenderType::line;
    else if (val == "tube")
      render_types[image_num][feature_num] = RenderType::tube;
    else
    {
      std::ostringstream message;
      message << "Invalid render type (" << val << ") in input file.";
      throw BlacklightException(message.str().c_str());
    }
  }

  // Read quantity
  else if (key.size() >= 12 and key.compare(key.size() - 9, key.npos, "_quantity") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    if (val == "rho")
      render_quantities[image_num][feature_num] = CellValues::rho;
    else if (val == "n_e")
      render_quantities[image_num][feature_num] = CellValues::n_e;
    else if (val == "p_gas")
      render_quantities[image_num][feature_num] = CellValues::p_gas;
    else if (val == "Theta_e")
      render_quantities[image_num][feature_num] = CellValues::theta_e;
    else if (val == "B")
      render_quantities[image_num][feature_num] = CellValues::bb;
    else if (val == "sigma")
      render_quantities[image_num][feature_num] = CellValues::sigma;
    else if (val == "beta_inverse")
      render_quantities[image_num][feature_num] = CellValues::beta_inv;
    else
    {
      std::ostringstream message;
      message << "Invalid render quantity (" << val << ") in input file.";
      throw BlacklightException(message.str().c_str());
    }
  }

  // Read minimum
  else if (key.size() >= 7 and key.compare(key.size() - 4, key.npos,  "_min") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_min_vals[image_num][feature_num] = std::stod(val);
  }

  // Read maximum
  else if (key.size() >= 7 and key.compare(key.size() - 4, key.npos, "_max") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_max_vals[image_num][feature_num] = std::stod(val);
  }

  // Read threshold
  else if (key.size() >= 10 and key.compare(key.size() - 7, key.npos,  "_thresh") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_thresh_vals[image_num][feature_num] = std::stod(val);
  }

  // Read radius
  else if (key.size() >= 5 and key.compare(key.size() - 2, key.npos,  "_r") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_r_vals[image_num][feature_num] = std::stod(val);
  }

  // Read tau scale
  else if (key.size() >= 13 and key.compare(key.size() - 10, key.npos, "_tau_scale") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_tau_scales[image_num][feature_num] = std::stod(val);
  }

  // Read opacity
  else if (key.size() >= 11 and key.compare(key.size() - 8, key.npos, "_opacity") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_opacities[image_num][feature_num] = std::stod(val);
  }

  // Read RGB values
  else if (key.size() >= 7 and key.compare(key.size() - 4, key.npos, "_rgb") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    double r, g, b;
    ReadTriple(val, &r, &g, &b);
    double x, y, z;
    RGBToXYZ(r, g, b, &x, &y, &z);
    render_x_vals[image_num][feature_num] = x;
    render_y_vals[image_num][feature_num] = y;
    render_z_vals[image_num][feature_num] = z;
  }

  // Read XYZ values
  else if (key.size() >= 7 and key.compare(key.size() - 4, key.npos, "_xyz") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    ReadTriple(val, &render_x_vals[image_num][feature_num], &render_y_vals[image_num][feature_num],
        &render_z_vals[image_num][feature_num]);
  }

  // Read light values
  else if (key.size() >= 9 and key.compare(key.size() - 6, key.npos, "_light") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    ReadTriple(val, &render_lx_vals[image_num][feature_num],
        &render_ly_vals[image_num][feature_num], &render_lz_vals[image_num][feature_num]);
  }

  // Read streamline filename
  else if (key.size() >= 8 and key.compare(key.size() - 5, key.npos, "_file") == 0)
  {
    std::size_t pos;
    int image_num = std::stoi(key, &pos) - 1;
    int feature_num = std::stoi(key.substr(pos + 1)) - 1;
    if (image_num >= render_num_images.value()
        or feature_num >= render_num_features[image_num].value())
      return;
    render_stream_files[image_num][feature_num] = val;
  }

  // Handle unknown entry
  else
  {
    std::ostringstream message;
    message << "Unknown key (render_" << key << ") in input file.";
    throw BlacklightException(message.str().c_str());
  }
  return;
}
