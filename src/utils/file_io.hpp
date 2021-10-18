// Blacklight file I/O header

#ifndef FILE_IO_H_
#define FILE_IO_H_

// C++ headers
#include <fstream>  // ifstream, ofstream

// Blacklight headers
#include "array.hpp"  // Array

//--------------------------------------------------------------------------------------------------

// Functions for writing binary data
void WriteBinary(std::ofstream *p_stream, int val);
void WriteBinary(std::ofstream *p_stream, double val);
void WriteBinary(std::ofstream *p_stream, double vals[], int num);
template<typename type> void WriteBinary(std::ofstream *p_stream, const Array<type> &array);

// Functions for reading binary data
void ReadBinary(std::ifstream *p_stream, int *p_val);
void ReadBinary(std::ifstream *p_stream, double *p_val);
void ReadBinary(std::ifstream *p_stream, double vals[], int num);
template<typename type> void ReadBinary(std::ifstream *p_stream, Array<type> *p_array);

#endif
