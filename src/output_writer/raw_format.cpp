// Blacklight output writer - raw output format

// Blacklight headers
#include "output_writer.hpp"
#include "../utils/array.hpp"    // Array
#include "../utils/file_io.hpp"  // WriteBinary

//--------------------------------------------------------------------------------------------------

// Output writer for a raw binary dump without metadata
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Writes data to file.
//   Output file endianness will match that of machine writing file.
void OutputWriter::WriteRaw()
{
  WriteBinary(p_output_stream, image[0].data, image[0].n_tot);
  return;
}
