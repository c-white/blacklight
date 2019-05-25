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

  // Data
  const std::string input_file;
  std::string data_file;
  double m;
  double a;

  // Functions
  void read();
  static bool removeable_space(unsigned char c);
};

#endif
