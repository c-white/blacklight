// Blacklight output writer - raw output format

// C++ headers
#include <fstream>  // ofstream
#include <ios>      // streamsize

// Blacklight headers
#include "output_writer.hpp"
#include "../utils/array.hpp"  // Array

//--------------------------------------------------------------------------------------------------

// Output writer for a raw binary dump without metadata
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
//   Output file endianness will match that of machine writing file.
void OutputWriter::WriteRaw()
{
  const char *data_pointer = reinterpret_cast<const char *>(image.data);
  std::streamsize data_size = static_cast<std::streamsize>(image.GetNumBytes());
  p_output_stream->write(data_pointer, data_size);
  return;
}
