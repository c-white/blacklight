// Blacklight multidimensional array

// C++ headers
#include <complex>  // complex
#include <cstddef>  // size_t
#include <cstring>  // memcpy
#include <limits>   // numeric_limits

// Blacklight headers
#include "array.hpp"
#include "exceptions.hpp"  // BlacklightException

// Instantiations
template struct Array<bool>;
template struct Array<int>;
template struct Array<float>;
template struct Array<double>;
template struct Array<std::complex<double>>;

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (empty)
// Inputs: (none)
template<typename type> Array<type>::Array() {}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (1D)
// Inputs:
//   n1_: size of only dimension
template<typename type> Array<type>::Array(int n1_)
  : n1(n1_),
    n2(1),
    n3(1),
    n4(1),
    n5(1)
{
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (2D)
// Inputs:
//   n2_: size of outermost dimension
//   n1_: size of innermost dimension
template<typename type> Array<type>::Array(int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(1),
    n4(1),
    n5(1)
{
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (3D)
// Inputs:
//   n3_: size of outermost dimension
//   n2_: size of intermediate dimension
//   n1_: size of innermost dimension
template<typename type> Array<type>::Array(int n3_, int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(n3_),
    n4(1),
    n5(1)
{
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (4D)
// Inputs:
//   n4_: size of outermost dimension
//   n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
template<typename type> Array<type>::Array(int n4_, int n3_, int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(n3_),
    n4(n4_),
    n5(1)
{
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array constructor (5D)
// Inputs:
//   n5_: size of outermost dimension
//   n4_, n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
template<typename type> Array<type>::Array(int n5_, int n4_, int n3_, int n2_, int n1_)
  : n1(n1_),
    n2(n2_),
    n3(n3_),
    n4(n4_),
    n5(n5_)
{
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array copy constructor
template<typename type> Array<type>::Array(const Array<type> &source)
{
  data = source.data;
  n1 = source.n1;
  n2 = source.n2;
  n3 = source.n3;
  n4 = source.n4;
  n5 = source.n5;
  n_tot = source.n_tot;
  allocated = source.allocated;
  is_copy = true;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array copy assignment constructor
template<typename type> Array<type> &Array<type>::operator=(const Array<type> &source)
{
  data = source.data;
  n1 = source.n1;
  n2 = source.n2;
  n3 = source.n3;
  n4 = source.n4;
  n5 = source.n5;
  n_tot = source.n_tot;
  allocated = source.allocated;
  is_copy = true;
  return *this;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array destructor
template<typename type> Array<type>::~Array()
{
  Deallocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (general)
// Inputs: (none)
// Outputs: (none)
template<typename type> void Array<type>::Allocate()
{
  if (allocated)
    throw BlacklightException("Attempting to reallocate array.");
  allocated = true;
  n_tot = n1 * n2 * n3 * n4 * n5;
  if (n_tot <= 0)
    throw BlacklightException("Attempting to allocate empty array.");
  data = new type[n_tot];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (1D)
// Inputs:
//   n1_: size of only dimension
// Outputs: (none)
template<typename type> void Array<type>::Allocate(int n1_)
{
  n1 = n1_;
  n2 = 1;
  n3 = 1;
  n4 = 1;
  n5 = 1;
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (2D)
// Inputs:
//   n2_: size of outermost dimension
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type> void Array<type>::Allocate(int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = 1;
  n4 = 1;
  n5 = 1;
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (3D)
// Inputs:
//   n3_: size of outermost dimension
//   n2_: size of intermediate dimension
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type> void Array<type>::Allocate(int n3_, int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  n4 = 1;
  n5 = 1;
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (4D)
// Inputs:
//   n4_: size of outermost dimension
//   n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type> void Array<type>::Allocate(int n4_, int n3_, int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  n4 = n4_;
  n5 = 1;
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array allocator (5D)
// Inputs:
//   n5_: size of outermost dimension
//   n4_, n3_, n2_: sizes of intermediate dimensions
//   n1_: size of innermost dimension
// Outputs: (none)
template<typename type> void Array<type>::Allocate(int n5_, int n4_, int n3_, int n2_, int n1_)
{
  n1 = n1_;
  n2 = n2_;
  n3 = n3_;
  n4 = n4_;
  n5 = n5_;
  Allocate();
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array deallocator
// Inputs: (none)
// Outputs: (none)
template<typename type> void Array<type>::Deallocate()
{
  if (allocated and not is_copy)
    delete[] data;
  allocated = false;
  return;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (1D)
// Inputs:
//   i1: only index
// Outputs:
//   returned value: element
template<typename type> type Array<type>::operator()(int i1) const
{
  return data[i1];
}
template<typename type> type Array<type>::operator()(unsigned int i1) const
{
  return data[i1];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (2D)
// Inputs:
//   i2: outermost index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type Array<type>::operator()(int i2, int i1) const
{
  return data[i1 + n1 * i2];
}
template<typename type> type Array<type>::operator()(unsigned int i2, unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * i2];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (3D)
// Inputs:
//   i3: outermost index
//   i2: intermediate index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type Array<type>::operator()(int i3, int i2, int i1) const
{
  return data[i1 + n1 * (i2 + n2 * i3)];
}
template<typename type> type Array<type>::operator()(unsigned int i3, unsigned int i2,
    unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2) * i3)];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (4D)
// Inputs:
//   i4: outermost index
//   i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type Array<type>::operator()(int i4, int i3, int i2, int i1) const
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * i4))];
}
template<typename type> type Array<type>::operator()(unsigned int i4, unsigned int i3,
    unsigned int i2, unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * i4))];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read accessors (5D)
// Inputs:
//   i5: outermost index
//   i4, i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type Array<type>::operator()(int i5, int i4, int i3, int i2, int i1) const
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * i5)))];
}
template<typename type> type Array<type>::operator()(unsigned int i5, unsigned int i4,
    unsigned int i3, unsigned int i2, unsigned int i1) const
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * (i4 + static_cast<unsigned int>(n4) * i5)))];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (1D)
// Inputs:
//   i1: only index
// Outputs:
//   returned value: element
template<typename type> type &Array<type>::operator()(int i1)
{
  return data[i1];
}
template<typename type> type &Array<type>::operator()(unsigned int i1)
{
  return data[i1];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (2D)
// Inputs:
//   i2: outermost index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type &Array<type>::operator()(int i2, int i1)
{
  return data[i1 + n1 * i2];
}
template<typename type> type &Array<type>::operator()(unsigned int i2, unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * i2];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (3D)
// Inputs:
//   i3: outermost index
//   i2: intermediate index
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type &Array<type>::operator()(int i3, int i2, int i1)
{
  return data[i1 + n1 * (i2 + n2 * i3)];
}
template<typename type> type &Array<type>::operator()(unsigned int i3, unsigned int i2,
    unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2) * i3)];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (4D)
// Inputs:
//   i4: outermost index
//   i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type &Array<type>::operator()(int i4, int i3, int i2, int i1)
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * i4))];
}
template<typename type> type &Array<type>::operator()(unsigned int i4, unsigned int i3,
    unsigned int i2, unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * i4))];
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array read/write accessors (5D)
// Inputs:
//   i5: outermost index
//   i4, i3, i2: intermediate indices
//   i1: innermost index
// Outputs:
//   returned value: element
template<typename type> type &Array<type>::operator()(int i5, int i4, int i3, int i2, int i1)
{
  return data[i1 + n1 * (i2 + n2 * (i3 + n3 * (i4 + n4 * i5)))];
}
template<typename type> type &Array<type>::operator()(unsigned int i5, unsigned int i4,
    unsigned int i3, unsigned int i2, unsigned int i1)
{
  return data[i1 + static_cast<unsigned int>(n1) * (i2 + static_cast<unsigned int>(n2)
      * (i3 + static_cast<unsigned int>(n3) * (i4 + static_cast<unsigned int>(n4) * i5)))];
}

//--------------------------------------------------------------------------------------------------

// Swap function
// Inputs:
//   other: array with which this should be swapped
// Outputs: (none)
template<typename type> void Array<type>::Swap(Array<type> &other)
{
  // Only allow swapping of allocated arrays
  if (not allocated or not other.allocated)
    throw BlacklightException("Attempting to swap an unallocated array.");

  // Make temporary shallow copy of other
  bool other_is_copy = other.is_copy;
  Array<type> temp(other);

  // Copy data from this to other
  other.data = data;
  other.n1 = n1;
  other.n2 = n2;
  other.n3 = n3;
  other.n4 = n4;
  other.n5 = n5;
  other.n_tot = n_tot;
  other.is_copy = is_copy;

  // Copy data from tmp to this
  data = temp.data;
  n1 = temp.n1;
  n2 = temp.n2;
  n3 = temp.n3;
  n4 = temp.n4;
  n5 = temp.n5;
  n_tot = temp.n_tot;
  is_copy = other_is_copy;
  return;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array slicing
// Inputs:
//   dimension: dimension of slice starting at 1 for innermost dimension
//   index: index of slice starting at 0 for first slice in given dimension
// Outputs: (none)
// Notes:
//   Moves data pointer and redefines data shape.
//   Checks to make sure only a shallow copy is being sliced so that destructor will not try to
//       deallocate memory.
template<typename type> void Array<type>::Slice(int dimension, int index)
{
  // Check that this is a shallow copy
  if (not is_copy)
    throw BlacklightException("Attempting to slice array that is not shallow copy.");

  // Take 1D slice
  if (dimension == 1)
  {
    if (index < 0 or index >= n1)
      throw BlacklightException("Attempting to slice outside array bounds.");
    data += index;
    n1 = 1;
    n2 = 1;
    n3 = 1;
    n4 = 1;
    n5 = 1;
  }

  // Take 2D slice
  else if (dimension == 2)
  {
    if (index < 0 or index >= n2)
      throw BlacklightException("Attempting to slice outside array bounds.");
    data += index * n1;
    n2 = 1;
    n3 = 1;
    n4 = 1;
    n5 = 1;
  }

  // Take 3D slice
  else if (dimension == 3)
  {
    if (index < 0 or index >= n3)
      throw BlacklightException("Attempting to slice outside array bounds.");
    data += index * n2 * n1;
    n3 = 1;
    n4 = 1;
    n5 = 1;
  }

  // Take 4D slice
  else if (dimension == 4)
  {
    if (index < 0 or index >= n4)
      throw BlacklightException("Attempting to slice outside array bounds.");
    data += index * n3 * n2 * n1;
    n4 = 1;
    n5 = 1;
  }

  // Take 5D slice
  else if (dimension == 5)
  {
    if (index < 0 or index >= n5)
      throw BlacklightException("Attempting to slice outside array bounds.");
    data += index * n4 * n3 * n2 * n1;
    n5 = 1;
  }

  // Account for invalid dimension
  else
    throw BlacklightException("Attempting to slice at invalid dimension.");

  // Recalculate number of elements
  n_tot = n1 * n2 * n3 * n4 * n5;
  return;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional array zeroing
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets any allocated array to be zero.
template<typename type> void Array<type>::Zero()
{
  if (allocated)
    for (int n = 0; n < n_tot; n++)
      data[n] = static_cast<type>(0);
  return;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional float array NaN-initializing
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets any allocated array to be NaN.
template<> void Array<float>::SetNaN()
{
  if (allocated)
    for (int n = 0; n < n_tot; n++)
      data[n] = std::numeric_limits<float>::quiet_NaN();
  return;
}

//--------------------------------------------------------------------------------------------------

// Multidimensional double array NaN-initializing
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Sets any allocated array to be NaN.
template<> void Array<double>::SetNaN()
{
  if (allocated)
    for (int n = 0; n < n_tot; n++)
      data[n] = std::numeric_limits<double>::quiet_NaN();
  return;
}

//--------------------------------------------------------------------------------------------------

// Byte size of array
// Inputs: (none)
// Outputs:
//   returned value: number of bytes in allocated array
template<typename type> std::size_t Array<type>::GetNumBytes() const
{
  return static_cast<std::size_t>(n_tot) * sizeof(type);
}
