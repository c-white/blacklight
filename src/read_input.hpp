// Ray Trace input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// C++ headers
#include <string>  // string

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

  // Data - image
  double im_r;
  double im_th;
  double im_ph;
  double im_rot;
  double im_width;
  int im_res;

  // Data - rays
  double ray_step;
  int ray_max_steps;
  bool flat;

  // Functions
  void Read();
  static bool RemoveableSpace(unsigned char c);
  bool ReadBool(const std::string &string);
};

#endif
