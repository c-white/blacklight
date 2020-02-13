// Blacklight input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// C++ headers
#include <string>  // string

// Blacklight headers
#include "blacklight.hpp"  // Coordinates

//--------------------------------------------------------------------------------------------------

// Input reader
struct InputReader
{
  // Constructor
  InputReader(const std::string input_file_);

  // Data - general
  const std::string input_file;

  // Data - file names
  std::string data_file;
  std::string output_file;

  // Data - coordinates
  double bh_m;
  double bh_a;
  Coordinates coord;

  // Data - units
  double m_msun;
  double rho_unit;

  // Data - plasma
  double plasma_mu;
  double plasma_ne_ni;
  double plasma_rat_high;
  double plasma_rat_low;
  double plasma_sigma_max;

  // Data - image
  double im_r;
  double im_th;
  double im_ph;
  double im_rot;
  double im_width;
  int im_res;
  double im_freq;
  bool im_pole;

  // Data - rays
  double ray_step;
  int ray_max_steps;
  bool ray_flat;

  // Data - performance
  int num_threads;

  // Functions
  void Read();
  static bool RemoveableSpace(unsigned char c);
  double ReadPole(const std::string &string, bool *p_pole_flag);
  bool ReadBool(const std::string &string);
  Coordinates ReadCoordinates(const std::string &string);
};

#endif
