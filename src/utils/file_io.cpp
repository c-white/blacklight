// Blacklight file I/O

// C++ headers
#include <cstddef>  // size_t
#include <fstream>  // ifstream, ofstream
#include <ios>      // streamsize

// Blacklight headers
#include "file_io.hpp"
#include "array.hpp"    // Array

// Instantiations
template void WriteBinary<int>(std::ofstream *p_stream, int val);
template void WriteBinary<double>(std::ofstream *p_stream, double vals[], long int num);
template void WriteBinary<bool>(std::ofstream *p_stream, const Array<bool> &array);
template void WriteBinary<int>(std::ofstream *p_stream, const Array<int> &array);
template void WriteBinary<double>(std::ofstream *p_stream, const Array<double> &array);
template void ReadBinary<int>(std::ifstream *p_stream, int *p_val);
template void ReadBinary<float>(std::ifstream *p_stream, float vals[], long int num);
template void ReadBinary<double>(std::ifstream *p_stream, double vals[], long int num);
template void ReadBinary<bool>(std::ifstream *p_stream, Array<bool> *p_array);
template void ReadBinary<int>(std::ifstream *p_stream, Array<int> *p_array);
template void ReadBinary<double>(std::ifstream *p_stream, Array<double> *p_array);

//--------------------------------------------------------------------------------------------------

// Function for writing single value to file
// Inputs:
//   *p_stream: open ofstream for file being written
//   val: integer to write
// Outputs: (none)
template<typename type> void WriteBinary(std::ofstream *p_stream, type val)
{
  const char *data_pointer = reinterpret_cast<const char *>(&val);
  std::streamsize data_size = sizeof(type);
  p_stream->write(data_pointer, data_size);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for writing array to file
// Inputs:
//   *p_stream: open ofstream for file being written
//   vals: array of doubles to write
//   num: number of elements
// Outputs: (none)
template<typename type> void WriteBinary(std::ofstream *p_stream, type vals[], long int num)
{
  const char *data_pointer = reinterpret_cast<const char *>(vals);
  std::streamsize data_size =
      static_cast<std::streamsize>(static_cast<std::size_t>(num) * sizeof(type));
  p_stream->write(data_pointer, data_size);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for writing a double Array to file
// Inputs:
//   *p_stream: open ofstream for file being written
//   array: Array of doubles to write
// Outputs: (none)
template<typename type> void WriteBinary(std::ofstream *p_stream, const Array<type> &array)
{
  WriteBinary(p_stream, array.n1);
  WriteBinary(p_stream, array.n2);
  WriteBinary(p_stream, array.n3);
  WriteBinary(p_stream, array.n4);
  WriteBinary(p_stream, array.n5);
  const char *data_pointer = reinterpret_cast<const char *>(array.data);
  std::streamsize data_size = static_cast<std::streamsize>(array.GetNumBytes());
  p_stream->write(data_pointer, data_size);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for reading single value from file
// Inputs:
//   *p_stream: open ifstream for file being read
// Outputs:
//   *p_val: value read from file
template<typename type> void ReadBinary(std::ifstream *p_stream, type *p_val)
{
  char *data_pointer = reinterpret_cast<char *>(p_val);
  std::streamsize data_size = sizeof(type);
  p_stream->read(data_pointer, data_size);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for reading array from file
// Inputs:
//   *p_stream: open ifstream for file being read
//   num: number of elements
// Outputs:
//   vals: array read from file
template<typename type> void ReadBinary(std::ifstream *p_stream, type vals[], long int num)
{
  char *data_pointer = reinterpret_cast<char *>(vals);
  std::streamsize data_size =
      static_cast<std::streamsize>(static_cast<std::size_t>(num) * sizeof(type));
  p_stream->read(data_pointer, data_size);
  return;
}

//--------------------------------------------------------------------------------------------------

// Function for reading Array from file
// Inputs:
//   *p_stream: open ifstream for file being read
// Outputs:
//   *p_array: Array of doubles read from file
template<typename type> void ReadBinary(std::ifstream *p_stream, Array<type> *p_array)
{
  int n1, n2, n3, n4, n5;
  ReadBinary(p_stream, &n1);
  ReadBinary(p_stream, &n2);
  ReadBinary(p_stream, &n3);
  ReadBinary(p_stream, &n4);
  ReadBinary(p_stream, &n5);
  p_array->Allocate(n5, n4, n3, n2, n1);
  char *data_pointer = reinterpret_cast<char *>(p_array->data);
  std::streamsize data_size = static_cast<std::streamsize>(p_array->GetNumBytes());
  p_stream->read(data_pointer, data_size);
  return;
}
