// Blacklight file I/O header

#ifndef FILE_IO_H_
#define FILE_IO_H_

// C++ headers
#include <fstream>  // ifstream, ofstream

// Blacklight headers
#include "array.hpp"  // Array

//--------------------------------------------------------------------------------------------------

// Functions for writing binary data
template<typename type> void WriteBinary(std::ofstream *p_stream, type val);
template<typename type> void WriteBinary(std::ofstream *p_stream, type vals[], long int num);
template<typename type> void WriteBinary(std::ofstream *p_stream, const Array<type> &array);

// Functions for reading binary data
template<typename type> void ReadBinary(std::ifstream *p_stream, type *p_val);
template<typename type> void ReadBinary(std::ifstream *p_stream, type vals[], long int num);
template<typename type> void ReadBinary(std::ifstream *p_stream, Array<type> *p_array);

#endif
