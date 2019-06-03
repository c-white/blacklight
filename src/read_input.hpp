// Ray Trace input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// C++ headers
#include <string>  // string

//--------------------------------------------------------------------------------------------------

// Input reader
struct input_reader
{
  // Constructor
  input_reader(const std::string input_file_);

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
  double im_step;

  // Functions
  void read();
  static bool removeable_space(unsigned char c);
};

#endif
