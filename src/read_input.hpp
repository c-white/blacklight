// Ray Trace input reader header

#ifndef READ_INPUT_H_
#define READ_INPUT_H_

// Input reader
struct input_reader
{
  // Constructor and destructor
  input_reader(const char *input_file_);
  ~input_reader();

  // Data
  const char *input_file;
  char *data_file;
  double m;
  double a;
};

#endif
