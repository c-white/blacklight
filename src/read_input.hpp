// Ray Trace input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// C++ headers
#include <string>  // string

// Input reader
struct input_reader
{
  // Constructor and destructor
  input_reader(const std::string input_file_);
  ~input_reader();

  // Data
  const std::string input_file;
  char *data_file;
  double m;
  double a;

  // Functions
  static bool removeable_space(unsigned char c);
};

#endif
