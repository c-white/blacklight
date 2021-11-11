// Blacklight output writer - NumPy output formats

// C++ headers
#include <cstdint>  // uint8_t, uint16_t
#include <cstdio>   // size_t, snprintf
#include <cstring>  // memcpy, memset
#include <fstream>  // ofstream
#include <ios>      // streamsize

// Blacklight headers
#include "output_writer.hpp"
#include "../blacklight.hpp"        // enums
#include "../utils/array.hpp"       // Array
#include "../utils/exceptions.hpp"  // BlacklightException

//--------------------------------------------------------------------------------------------------

// Output writer for a NumPy .npy file
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
void OutputWriter::WriteNpy()
{
  uint8_t *image_buffer;
  std::size_t image_size = GenerateNpyFromArray(image[0], 3, &image_buffer);
  char *buffer = reinterpret_cast<char *>(image_buffer);
  std::streamsize buffer_length = static_cast<std::streamsize>(image_size);
  p_output_stream->write(buffer, buffer_length);
  delete[] image_buffer;
  return;
}

//--------------------------------------------------------------------------------------------------

// Output writer for a NumPy .npz file
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
//   NumPy .npz files are simply ZIP files with each record containing a .npy file and named
//       according to the name of the array.
//   Implements ZIP format version 2.0: pkware.cachefly.net/webdocs/casestudies/APPNOTE.TXT
//   Does not implement ZIP64, and thus has a reasonably low limit to data sizes that should still
//       be plenty for almost all ray tracing applications.
void OutputWriter::WriteNpz()
{
  // Prepare buffers for data and headers
  int num_image_arrays = (image_light ? 1 : 0)
      + (image_light and model_type == ModelType::simulation and image_polarization ? 3 : 0)
      + (image_time ? 1 : 0) + (image_length ? 1 : 0) + (image_lambda ? 1 : 0)
      + (image_emission ? 1 : 0) + (image_tau ? 1 : 0)
      + (image_lambda_ave ? CellValues::num_cell_values : 0)
      + (image_emission_ave ? CellValues::num_cell_values : 0)
      + (image_tau_int ? CellValues::num_cell_values : 0);
  int num_full_arrays =
      (output_camera ? 1 : 0) + num_image_arrays + (render_num_images > 0 ? 1 : 0);
  int num_arrays = 3 + 1 + (adaptive_max_level > 0 ? 1 : 0) + num_full_arrays
       + adaptive_num_levels_array(0) * (1 + num_full_arrays);
  const int max_name_length = 128;
  char *name_buffer = new char[max_name_length];
  uint8_t **data_buffers = new uint8_t *[num_arrays];
  uint8_t **local_header_buffers = new uint8_t *[num_arrays];
  uint8_t **central_header_buffers = new uint8_t *[num_arrays];
  std::size_t *data_lengths = new std::size_t[num_arrays];
  std::size_t *local_header_lengths = new std::size_t[num_arrays];
  std::size_t *central_header_lengths = new std::size_t[num_arrays];
  int array_offset = 0;

  // Write output parameters and metadata to buffers
  data_lengths[array_offset] =
      GenerateNpyFromArray(mass_msun_array, 1, &data_buffers[array_offset]);
  local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
      data_lengths[array_offset], "mass_msun", &local_header_buffers[array_offset]);
  array_offset++;
  data_lengths[array_offset] =
      GenerateNpyFromArray(camera_width_array, 1, &data_buffers[array_offset]);
  local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
      data_lengths[array_offset], "width", &local_header_buffers[array_offset]);
  array_offset++;
  data_lengths[array_offset] =
      GenerateNpyFromArray(image_frequencies, 1, &data_buffers[array_offset]);
  local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
      data_lengths[array_offset], "frequency", &local_header_buffers[array_offset]);
  array_offset++;

  // Write number of adaptive levels and metadata to buffers
  data_lengths[array_offset] =
      GenerateNpyFromArray(adaptive_num_levels_array, 1, &data_buffers[array_offset]);
  local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
      data_lengths[array_offset], "adaptive_num_levels", &local_header_buffers[array_offset]);
  array_offset++;

  // Write numbers of adaptive blocks and metadata to buffers
  if (adaptive_max_level > 0)
  {
    data_lengths[array_offset] =
        GenerateNpyFromArray(block_counts_array, 1, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "adaptive_num_blocks", &local_header_buffers[array_offset]);
    array_offset++;
  }

  // Write root camera data and metadata to buffers
  if (output_camera)
  {
    if (camera_type == Camera::plane)
    {
      data_lengths[array_offset] =
          GenerateNpyFromArray(camera_pos[0], 3, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], "positions", &local_header_buffers[array_offset]);
    }
    else if (camera_type == Camera::pinhole)
    {
      data_lengths[array_offset] =
          GenerateNpyFromArray(camera_dir[0], 3, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], "directions", &local_header_buffers[array_offset]);
    }
    array_offset++;
  }

  // Write root intensity image data and metadata to buffers
  Array<double> image_deep_copy;
  if (image_light or image_lambda_ave or image_emission_ave or image_tau_int)
    image_deep_copy.Allocate(image_num_frequencies, camera_resolution, camera_resolution);
  int num_pix = camera_resolution * camera_resolution;
  int num_dims = image_num_frequencies == 1 ? 2 : 3;
  int image_stride =
      image_light and model_type == ModelType::simulation and image_polarization ? 4 : 1;
  if (image_light)
  {
    for (int l = 0; l < image_num_frequencies; l++)
      image_deep_copy.CopyFrom(image[0], l * image_stride * num_pix, l * num_pix, num_pix);
    data_lengths[array_offset] =
        GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "I_nu", &local_header_buffers[array_offset]);
    array_offset++;
    if (model_type == ModelType::simulation and image_polarization)
    {
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[0], (l * image_stride + 1) * num_pix, l * num_pix, num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], "Q_nu", &local_header_buffers[array_offset]);
      array_offset++;
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[0], (l * image_stride + 2) * num_pix, l * num_pix, num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], "U_nu", &local_header_buffers[array_offset]);
      array_offset++;
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[0], (l * image_stride + 3) * num_pix, l * num_pix, num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], "V_nu", &local_header_buffers[array_offset]);
      array_offset++;
    }
  }

  // Write root alternate image data and metadata to buffers
  Array<double> image_shallow_copy;
  if (image_time)
  {
    image_shallow_copy = image[0];
    image_shallow_copy.Slice(3, image_offset_time, image_offset_time);
    data_lengths[array_offset] =
        GenerateNpyFromArray(image_shallow_copy, 2, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "time", &local_header_buffers[array_offset]);
    array_offset++;
  }
  if (image_length)
  {
    image_shallow_copy = image[0];
    image_shallow_copy.Slice(3, image_offset_length, image_offset_length);
    data_lengths[array_offset] =
        GenerateNpyFromArray(image_shallow_copy, 2, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "length", &local_header_buffers[array_offset]);
    array_offset++;
  }
  if (image_lambda)
  {
    image_shallow_copy = image[0];
    image_shallow_copy.Slice(3, image_offset_lambda,
        image_offset_lambda + image_num_frequencies - 1);
    data_lengths[array_offset] =
        GenerateNpyFromArray(image_shallow_copy, num_dims, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "lambda", &local_header_buffers[array_offset]);
    array_offset++;
  }
  if (image_emission)
  {
    image_shallow_copy = image[0];
    image_shallow_copy.Slice(3, image_offset_emission,
        image_offset_emission + image_num_frequencies - 1);
    data_lengths[array_offset] =
        GenerateNpyFromArray(image_shallow_copy, num_dims, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "emission", &local_header_buffers[array_offset]);
    array_offset++;
  }
  if (image_tau)
  {
    image_shallow_copy = image[0];
    image_shallow_copy.Slice(3, image_offset_tau, image_offset_tau + image_num_frequencies - 1);
    data_lengths[array_offset] =
        GenerateNpyFromArray(image_shallow_copy, num_dims, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "tau", &local_header_buffers[array_offset]);
    array_offset++;
  }
  if (image_lambda_ave)
    for (int n = 0; n < CellValues::num_cell_values; n++)
    {
      int num_written = std::snprintf(name_buffer, max_name_length, "lambda_ave_%s", cell_names[n]);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[0],
            (image_offset_lambda_ave + l * CellValues::num_cell_values + n) * num_pix, l * num_pix,
            num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
  if (image_emission_ave)
    for (int n = 0; n < CellValues::num_cell_values; n++)
    {
      int num_written =
          std::snprintf(name_buffer, max_name_length, "emission_ave_%s", cell_names[n]);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[0],
            (image_offset_emission_ave + l * CellValues::num_cell_values + n) * num_pix,
            l * num_pix, num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
  if (image_tau_int)
    for (int n = 0; n < CellValues::num_cell_values; n++)
    {
      int num_written = std::snprintf(name_buffer, max_name_length, "tau_int_%s", cell_names[n]);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[0],
            (image_offset_tau_int + l * CellValues::num_cell_values + n) * num_pix, l * num_pix,
            num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
  image_deep_copy.Deallocate();

  // Write root render data and metadata to buffers
  if (render_num_images > 0)
  {
    data_lengths[array_offset] = GenerateNpyFromArray(render[0], 4, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], "rendering", &local_header_buffers[array_offset]);
    array_offset++;
  }

  // Write adaptive data and metadata to buffers
  num_dims++;
  for (int level = 1; level <= adaptive_num_levels_array(0); level++)
  {
    // Write adaptive refinement structure and metadata to buffers
    int num_written = std::snprintf(name_buffer, max_name_length, "adaptive_block_locs_%d", level);
    if (num_written < 0 or num_written >= max_name_length)
      throw BlacklightException("Error naming output array.");
    data_lengths[array_offset] =
        GenerateNpyFromArray(camera_loc[level], 2, &data_buffers[array_offset]);
    local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
        data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
    array_offset++;

    // Write adaptive camera data and metadata to buffers
    if (output_camera)
    {
      if (camera_type == Camera::plane)
      {
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_positions_%d", level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        data_lengths[array_offset] =
            GenerateNpyFromArray(camera_pos[level], 4, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      }
      else if (camera_type == Camera::pinhole)
      {
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_directions_%d", level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        data_lengths[array_offset] =
            GenerateNpyFromArray(camera_dir[level], 4, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      }
      array_offset++;
    }

    // Write adaptive intensity image data and metadata to buffers
    if (image_light or image_lambda_ave or image_emission_ave or image_tau_int)
      image_deep_copy.Allocate(image_num_frequencies, block_counts_array(level),
          adaptive_block_size, adaptive_block_size);
    num_pix = block_counts_array(level) * adaptive_block_size * adaptive_block_size;
    if (image_light)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_I_nu_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      for (int l = 0; l < image_num_frequencies; l++)
        image_deep_copy.CopyFrom(image[level], l * image_stride * num_pix, l * num_pix, num_pix);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
      if (model_type == ModelType::simulation and image_polarization)
      {
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_Q_nu_%d", level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        for (int l = 0; l < image_num_frequencies; l++)
          image_deep_copy.CopyFrom(image[level], (l * image_stride + 1) * num_pix, l * num_pix,
              num_pix);
        data_lengths[array_offset] =
            GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
        array_offset++;
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_U_nu_%d", level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        for (int l = 0; l < image_num_frequencies; l++)
          image_deep_copy.CopyFrom(image[level], (l * image_stride + 2) * num_pix, l * num_pix,
              num_pix);
        data_lengths[array_offset] =
            GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
        array_offset++;
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_V_nu_%d", level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        for (int l = 0; l < image_num_frequencies; l++)
          image_deep_copy.CopyFrom(image[level], (l * image_stride + 3) * num_pix, l * num_pix,
              num_pix);
        data_lengths[array_offset] =
            GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
        array_offset++;
      }
    }

    // Write adaptive alternate image data and metadata to buffers
    if (image_time)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_time_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      image_shallow_copy = image[level];
      image_shallow_copy.Slice(4, image_offset_time, image_offset_time);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_shallow_copy, 3, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
    if (image_length)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_length_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      image_shallow_copy = image[level];
      image_shallow_copy.Slice(4, image_offset_length, image_offset_length);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_shallow_copy, 3, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
    if (image_lambda)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_lambda_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      image_shallow_copy = image[level];
      image_shallow_copy.Slice(4, image_offset_lambda,
          image_offset_lambda + image_num_frequencies - 1);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_shallow_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
    if (image_emission)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_emission_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      image_shallow_copy = image[level];
      image_shallow_copy.Slice(4, image_offset_emission,
          image_offset_emission + image_num_frequencies - 1);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_shallow_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
    if (image_tau)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_tau_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      image_shallow_copy = image[level];
      image_shallow_copy.Slice(4, image_offset_tau, image_offset_tau + image_num_frequencies - 1);
      data_lengths[array_offset] =
          GenerateNpyFromArray(image_shallow_copy, num_dims, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
    if (image_lambda_ave)
      for (int n = 0; n < CellValues::num_cell_values; n++)
      {
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_lambda_ave_%s_%d",
            cell_names[n], level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        for (int l = 0; l < image_num_frequencies; l++)
          image_deep_copy.CopyFrom(image[level],
              (image_offset_lambda_ave + l * CellValues::num_cell_values + n) * num_pix,
              l * num_pix, num_pix);
        data_lengths[array_offset] =
            GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
        array_offset++;
      }
    if (image_emission_ave)
      for (int n = 0; n < CellValues::num_cell_values; n++)
      {
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_emission_ave_%s_%d",
            cell_names[n], level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        for (int l = 0; l < image_num_frequencies; l++)
          image_deep_copy.CopyFrom(image[level],
              (image_offset_emission_ave + l * CellValues::num_cell_values + n) * num_pix,
              l * num_pix, num_pix);
        data_lengths[array_offset] =
            GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
        array_offset++;
      }
    if (image_tau_int)
      for (int n = 0; n < CellValues::num_cell_values; n++)
      {
        num_written = std::snprintf(name_buffer, max_name_length, "adaptive_tau_int_%s_%d",
            cell_names[n], level);
        if (num_written < 0 or num_written >= max_name_length)
          throw BlacklightException("Error naming output array.");
        for (int l = 0; l < image_num_frequencies; l++)
          image_deep_copy.CopyFrom(image[level],
              (image_offset_tau_int + l * CellValues::num_cell_values + n) * num_pix, l * num_pix,
              num_pix);
        data_lengths[array_offset] =
            GenerateNpyFromArray(image_deep_copy, num_dims, &data_buffers[array_offset]);
        local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
            data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
        array_offset++;
      }
    image_deep_copy.Deallocate();

    // Write adaptive render data and metadata to buffers
    if (render_num_images > 0)
    {
      num_written = std::snprintf(name_buffer, max_name_length, "adaptive_rendering_%d", level);
      if (num_written < 0 or num_written >= max_name_length)
        throw BlacklightException("Error naming output array.");
      data_lengths[array_offset] =
          GenerateNpyFromArray(render[level], 5, &data_buffers[array_offset]);
      local_header_lengths[array_offset] = GenerateZIPLocalFileHeader(data_buffers[array_offset],
          data_lengths[array_offset], name_buffer, &local_header_buffers[array_offset]);
      array_offset++;
    }
  }

  // Write central directory headers to buffers
  std::size_t offset = 0;
  for (int n = 0; n < num_arrays; n++)
  {
    central_header_lengths[n] = GenerateZIPCentralDirectoryHeader(local_header_buffers[n], offset,
        &central_header_buffers[n]);
    offset += local_header_lengths[n] + data_lengths[n];
  }

  // Write end of central directory record to buffer
  std::size_t central_directory_length = 0;
  for (int n = 0; n < num_arrays; n++)
    central_directory_length += central_header_lengths[n];
  uint8_t *end_of_directory_buffer = nullptr;
  std::size_t end_of_directory_length = GenerateZIPEndOfCentralDirectoryRecord(offset,
      central_directory_length, num_arrays, &end_of_directory_buffer);

  // Write buffers to file
  for (int n = 0; n < num_arrays; n++)
  {
    char *buffer = reinterpret_cast<char *>(local_header_buffers[n]);
    std::streamsize length = static_cast<std::streamsize>(local_header_lengths[n]);
    p_output_stream->write(buffer, length);
    buffer = reinterpret_cast<char *>(data_buffers[n]);
    length = static_cast<std::streamsize>(data_lengths[n]);
    p_output_stream->write(buffer, length);
  }
  for (int n = 0; n < num_arrays; n++)
  {
    char *buffer = reinterpret_cast<char *>(central_header_buffers[n]);
    std::streamsize length = static_cast<std::streamsize>(central_header_lengths[n]);
    p_output_stream->write(buffer, length);
  }
  char *buffer = reinterpret_cast<char *>(end_of_directory_buffer);
  std::streamsize length = static_cast<std::streamsize>(end_of_directory_length);
  p_output_stream->write(buffer, length);

  // Free memory
  for (int n = 0; n < num_arrays; n++)
  {
    delete[] data_buffers[n];
    delete[] local_header_buffers[n];
    delete[] central_header_buffers[n];
  }
  delete[] name_buffer;
  delete[] data_buffers;
  delete[] local_header_buffers;
  delete[] central_header_buffers;
  delete[] end_of_directory_buffer;
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for populating a buffer with the contents of a NumPy .npy file from a double Array
// Inputs:
//   array: Array to be written to buffer
//   num_dims: number of dimensions in Array
// Outputs:
//   p_buffer: value set to newly allocated buffer that has been initialized with file contents
//   returned value: length of allocated buffer
// Notes:
//   Must be run on little-endian machine.
//   If num_dims is greater than the number of non-singleton dimensions in the Array, the output
//       array will have leading singleton dimensions.
//   Implements NumPy .npy format version 1.0:
//       numpy.org/doc/stable/reference/generated/numpy.lib.format.html#module-numpy.lib.format
//   A .npy file is simply a binary dump of a NumPy array prepended with an ASCII string that can be
//       interpreted as a Python dictionary literal with certain entries.
//   Also overloaded for int Array.
std::size_t OutputWriter::GenerateNpyFromArray(const Array<double> &array, int num_dims,
    uint8_t **p_buffer)
{
  // Prepare buffer
  const std::size_t header_length = 128;
  std::size_t data_length = array.GetNumBytes();
  std::size_t buffer_length = header_length + data_length;
  *p_buffer = new uint8_t[buffer_length];

  // Write magic string and version number to buffer
  std::size_t length = 0;
  const char *magic_string = "\x93NUMPY";
  std::memcpy(*p_buffer + length, magic_string, 6);
  length += 6;
  const char *version = "\x01\x00";
  std::memcpy(*p_buffer + length, version, 2);
  length += 2;

  // Write length of header proper to buffer
  uint16_t header_len = static_cast<uint16_t>(header_length - length - 2);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&header_len), 2);
  length += 2;

  // Write header proper to buffer
  char *buffer_address = reinterpret_cast<char *>(*p_buffer + length);
  int num_written = -1;
  if (num_dims == 1 and array.n2 == 1 and array.n3 == 1 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d,)}", array.n1);
  else if (num_dims == 2 and array.n3 == 1 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d)}", array.n2, array.n1);
  else if (num_dims == 3 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d, %d)}", array.n3, array.n2,
        array.n1);
  else if (num_dims == 4 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d, %d, %d)}", array.n4, array.n3,
        array.n2, array.n1);
  else if (num_dims == 5)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d, %d, %d, %d)}", array.n5,
        array.n4, array.n3, array.n2, array.n1);
  else
    throw BlacklightException("Attempt to truncate array while writing to .npy format.");
  if (num_written < 0 or num_written > static_cast<int>(header_length - length - 1))
    throw BlacklightException("Error converting data to .npy format.");
  std::size_t num_spaces =
      static_cast<std::size_t>(static_cast<int>(header_length - length - 1) - num_written);
  if (num_spaces > 0)
    std::memset(buffer_address + num_written, ' ', num_spaces);
  (*p_buffer)[header_length-1] = '\n';
  length = header_length;

  // Write array to buffer
  const uint8_t *data_pointer = reinterpret_cast<const uint8_t *>(array.data);
  std::memcpy(*p_buffer + length, data_pointer, data_length);
  length += data_length;
  return buffer_length;
}

//--------------------------------------------------------------------------------------------------

// Function for populating a buffer with the contents of a NumPy .npy file from an int Array
// Inputs:
//   array: Array to be written to buffer
//   num_dims: number of dimensions in Array
// Outputs:
//   p_buffer: value set to newly allocated buffer that has been initialized with file contents
//   returned value: length of allocated buffer
// Notes:
//   Must be run on little-endian machine.
//   If num_dims is greater than the number of non-singleton dimensions in the Array, the output
//       array will have leading singleton dimensions.
//   Implements NumPy .npy format version 1.0:
//       numpy.org/doc/stable/reference/generated/numpy.lib.format.html#module-numpy.lib.format
//   A .npy file is simply a binary dump of a NumPy array prepended with an ASCII string that can be
//       interpreted as a Python dictionary literal with certain entries.
//   Also overloaded for double Array.
std::size_t OutputWriter::GenerateNpyFromArray(const Array<int> &array, int num_dims,
    uint8_t **p_buffer)
{
  // Prepare buffer
  const std::size_t header_length = 128;
  std::size_t data_length = array.GetNumBytes();
  std::size_t buffer_length = header_length + data_length;
  *p_buffer = new uint8_t[buffer_length];

  // Write magic string and version number to buffer
  std::size_t length = 0;
  const char *magic_string = "\x93NUMPY";
  std::memcpy(*p_buffer + length, magic_string, 6);
  length += 6;
  const char *version = "\x01\x00";
  std::memcpy(*p_buffer + length, version, 2);
  length += 2;

  // Write length of header proper to buffer
  uint16_t header_len = static_cast<uint16_t>(header_length - length - 2);
  std::memcpy(*p_buffer + length, reinterpret_cast<const uint8_t *>(&header_len), 2);
  length += 2;

  // Write header proper to buffer
  char *buffer_address = reinterpret_cast<char *>(*p_buffer + length);
  int num_written = -1;
  if (num_dims == 1 and array.n2 == 1 and array.n3 == 1 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<i4', 'fortran_order': False, 'shape': (%d,)}", array.n1);
  else if (num_dims == 2 and array.n3 == 1 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<i4', 'fortran_order': False, 'shape': (%d, %d)}", array.n2, array.n1);
  else if (num_dims == 3 and array.n4 == 1 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<i4', 'fortran_order': False, 'shape': (%d, %d, %d)}", array.n3, array.n2,
        array.n1);
  else if (num_dims == 4 and array.n5 == 1)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<i4', 'fortran_order': False, 'shape': (%d, %d, %d, %d)}", array.n4, array.n3,
        array.n2, array.n1);
  else if (num_dims == 5)
    num_written = std::snprintf(buffer_address, header_length - length,
        "{'descr': '<i4', 'fortran_order': False, 'shape': (%d, %d, %d, %d, %d)}", array.n5,
        array.n4, array.n3, array.n2, array.n1);
  else
    throw BlacklightException("Attempt to truncate array while writing to .npy format.");
  if (num_written < 0 or num_written > static_cast<int>(header_length - length - 1))
    throw BlacklightException("Error converting data to .npy format.");
  std::size_t num_spaces =
      static_cast<std::size_t>(static_cast<int>(header_length - length - 1) - num_written);
  if (num_spaces > 0)
    std::memset(buffer_address + num_written, ' ', num_spaces);
  (*p_buffer)[header_length-1] = '\n';
  length = header_length;

  // Write array to buffer
  const uint8_t *data_pointer = reinterpret_cast<const uint8_t *>(array.data);
  std::memcpy(*p_buffer + length, data_pointer, data_length);
  length += data_length;
  return buffer_length;
}
